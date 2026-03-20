"""
ProteinClaw - 蛋白质注释整合工具
================================

把"一个蛋白的所有已知注释"系统化取回并结构化输出
"""

from flask import Flask, render_template, request, jsonify
import requests
import re
from typing import Optional, Dict, List, Any, Tuple

app = Flask(__name__)
app.config['SECRET_KEY'] = 'proteinclaw-secret-key'


class ProteinAnnotator:
    """蛋白质注释整合器"""

    def __init__(self):
        self.uniprot_base = "https://rest.uniprot.org/uniprotkb"
        self.mygeneinfo_base = "https://mygene.info/v3"
        self.string_base = "https://string-db.org/api"
        self.pfam_base = "https://pfam.xfam.org"

    def normalize_id(self, input_id: str, species: str = "human") -> Dict[str, Any]:
        """
        Step 0: ID 规范化与映射

        支持输入：
        - UniProt Accession (Q9H3U1)
        - Gene symbol (SMOC1)
        - Ensembl ID (ENSG00000123456)

        返回标准化的 ID 映射
        """
        input_id = input_id.strip()

        # 判断输入类型
        if re.match(r'^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$', input_id):
            input_type = "uniprot"
        elif re.match(r'^ENSG[0-9]{11}$', input_id):
            input_type = "ensembl"
        else:
            input_type = "gene_symbol"

        result = {
            "input": input_id,
            "input_type": input_type,
            "species": species,
            "uniprot_id": None,
            "gene_symbol": None,
            "ensembl_id": None,
            "success": False
        }

        # 从 UniProt 直接查询
        if input_type == "uniprot":
            return self._fetch_uniprot_mapping(input_id, result)

        # 使用 UniProt search 直接搜索基因名
        try:
            if species == "human":
                taxon_id = "9606"
            elif species == "mouse":
                taxon_id = "10090"
            else:
                taxon_id = "9606"

            # 简化查询：只搜索输入内容
            url = f"https://rest.uniprot.org/uniprotkb/search?query={input_id}&format=json&size=5"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                if "results" in data and len(data["results"]) > 0:
                    # 找到匹配物种的第一个结果
                    for first_result in data["results"]:
                        if str(first_result.get("taxonomy", {}).get("taxonId", "")) == taxon_id:
                            uniprot_id = first_result["primaryAccession"]
                            # 找到 UniProt 后再获取详细映射
                            return self._fetch_uniprot_mapping(uniprot_id, result)
                    # 如果没找到对应物种，取第一个结果
                    first_result = data["results"][0]
                    uniprot_id = first_result["primaryAccession"]
                    return self._fetch_uniprot_mapping(uniprot_id, result)

        except Exception as e:
            print(f"ID mapping search error: {e}")

        return result

    def _fetch_uniprot_mapping(self, uniprot_id: str, result: Dict) -> Dict[str, Any]:
        """从 UniProt 获取映射信息"""
        try:
            url = f"{self.uniprot_base}/{uniprot_id}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                result["uniprot_id"] = uniprot_id

                # 提取基因名
                if "genes" in data and len(data["genes"]) > 0:
                    result["gene_symbol"] = data["genes"][0]["geneName"]

                # 提取 Ensembl ID
                if "dbReferences" in data:
                    for ref in data["dbReferences"]:
                        if ref["type"] == "Ensembl":
                            result["ensembl_id"] = ref["id"]
                            break

                result["success"] = True
        except Exception as e:
            print(f"UniProt mapping error: {e}")

        return result

    def fetch_basic_info(self, uniprot_id: str) -> Dict[str, Any]:
        """
        Step 1: 获取基础信息（UniProt 核心字段）
        """
        result = {
            "source": "UniProt",
            "protein_name": None,
            "gene_name": None,
            "length": None,
            "sequence": None,
            "reviewed": False,
            "subcellular_location": [],
            "function": None,
            "keywords": []
        }

        try:
            url = f"{self.uniprot_base}/{uniprot_id}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()

                result["protein_name"] = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
                if "genes" in data and len(data["genes"]) > 0:
                    result["gene_name"] = data["genes"][0]["geneName"]
                result["length"] = data.get("sequence", {}).get("length")
                result["sequence"] = data.get("sequence", {}).get("value")
                result["reviewed"] = data.get("reviewed", False)

                # 亚细胞定位
                if "comments" in data:
                    for comment in data["comments"]:
                        if comment["commentType"] == "SUBCELLULAR LOCATION":
                            if "locations" in comment:
                                for loc in comment["locations"]:
                                    if "location" in loc:
                                        result["subcellular_location"].append({
                                            "value": loc["location"]["value"],
                                            "source": "UniProt"
                                        })

                # 功能描述
                if "comments" in data:
                    for comment in data["comments"]:
                        if comment["commentType"] == "FUNCTION":
                            texts = comment.get("texts", [])
                            if texts:
                                result["function"] = texts[0]["value"]
                            break

                # 关键词
                if "keywords" in data:
                    for kw in data["keywords"]:
                        result["keywords"].append(kw["value"])

        except Exception as e:
            print(f"Basic info error: {e}")

        return result

    def fetch_go_annotations(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Step 2: 获取 GO 注释
        """
        result = []

        try:
            url = f"{self.uniprot_base}/{uniprot_id}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()

                if "uniProtKBCrossReferences" in data:
                    for ref in data["uniProtKBCrossReferences"]:
                        if ref["database"] == "GO":
                            go_id = ref["id"]
                            properties = {p["key"]: p["value"] for p in ref["properties"]}
                            category = properties.get("GoTerm")
                            evidence = properties.get("EvidenceCode", "IEA")
                            term_name = properties.get("GoName", "")

                            if category == "P":
                                category = "BP"
                            elif category == "F":
                                category = "MF"
                            elif category == "C":
                                category = "CC"

                            result.append({
                                "id": go_id,
                                "name": term_name,
                                "category": category,
                                "evidence": evidence,
                                "source": "Gene Ontology"
                            })

        except Exception as e:
            print(f"GO annotation error: {e}")

        return result

    def fetch_pathways(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Step 3: 获取 Pathway 注释（KEGG + Reactome）
        """
        result = []

        try:
            url = f"{self.uniprot_base}/{uniprot_id}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()

                if "uniProtKBCrossReferences" in data:
                    for ref in data["uniProtKBCrossReferences"]:
                        if ref["database"] == "KEGG":
                            result.append({
                                "id": ref["id"],
                                "name": ref.get("properties", [{}])[0].get("value", "") if ref.get("properties") else "",
                                "source": "KEGG"
                            })
                        elif ref["database"] == "Reactome":
                            result.append({
                                "id": ref["id"],
                                "name": ref.get("properties", [{}])[0].get("value", "") if ref.get("properties") else "",
                                "source": "Reactome"
                            })

        except Exception as e:
            print(f"Pathway error: {e}")

        return result

    def fetch_interactions_string(self, uniprot_id: str, species: str = "human", limit: int = 50) -> List[Dict[str, Any]]:
        """
        Step 4: 获取蛋白质相互作用（STRING）
        """
        result = []

        if species == "human":
            species_id = 9606
        elif species == "mouse":
            species_id = 10090
        else:
            species_id = 9606

        try:
            # Step 1: resolve ID to get preferred name
            url = f"{self.string_base}/json/resolve?identifiers={uniprot_id}&species={species_id}"
            resp = requests.get(url, timeout=10)
            if resp.status_code != 200:
                return result

            resolve_data = resp.json()
            if len(resolve_data) == 0:
                return result

            # Get the matched preferred name
            query_name = resolve_data[0].get("preferredName", uniprot_id)

            # Step 2: get network
            url = f"{self.string_base}/json/network?identifiers={uniprot_id}&species={species_id}&limit={limit}"
            resp = requests.get(url, timeout=15)
            if resp.status_code == 200:
                interactions = resp.json()

                for inter in interactions:
                    # STRING 返回的是两个蛋白
                    protein1 = inter.get("preferredName_A", "")
                    protein2 = inter.get("preferredName_B", "")

                    if protein1 == query_name:
                        partner_name = protein2
                    elif protein2 == query_name:
                        partner_name = protein1
                    else:
                        # 如果名字不匹配，用另一种方式：只要包含我们的蛋白就算
                        continue

                    score = inter.get("score", 0)
                    # 判断证据类型 - experimental score is the evidence from experiments
                    has_experimental = inter.get("escore", 0) > 0.1
                    evidence_type = "experimental" if has_experimental else "predicted"

                    result.append({
                        "partner": partner_name,
                        "score": score,
                        "evidence": evidence_type,
                        "source": "STRING"
                    })

                # 按分数排序，去重
                seen = set()
                result = [x for x in result if not (x["partner"] in seen or seen.add(x["partner"]))]
                result.sort(key=lambda x: x["score"], reverse=True)

                # 限制数量
                if len(result) > limit:
                    result = result[:limit]

        except Exception as e:
            print(f"STRING interaction error: {e}")

        return result

    def fetch_domains(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Step 5: 获取结构域注释（Pfam + InterPro）
        """
        result = []

        try:
            url = f"{self.uniprot_base}/{uniprot_id}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()

                if "features" in data:
                    for feature in data["features"]:
                        if feature["type"] in ["Domain", "Pfam", "InterPro"]:
                            source = feature.get("type", "Domain")
                            if source == "Domain":
                                source = "UniProt"

                            result.append({
                                "name": feature.get("description", ""),
                                "start": feature.get("location", {}).get("start", {}).get("value"),
                                "end": feature.get("location", {}).get("end", {}).get("value"),
                                "source": source
                            })

                # 从 crossReferences 获取 Pfam
                if "uniProtKBCrossReferences" in data:
                    for ref in data["uniProtKBCrossReferences"]:
                        if ref["database"] == "Pfam":
                            # 位置信息在 UniProt features 里已经有了
                            result.append({
                                "name": ref["id"] + " - " + (ref.get("properties", [{}])[0].get("value", "") if ref.get("properties") else ""),
                                "start": None,
                                "end": None,
                                "source": "Pfam"
                            })
                        elif ref["database"] == "InterPro":
                            result.append({
                                "name": ref["id"] + " - " + (ref.get("properties", [{}])[0].get("value", "") if ref.get("properties") else ""),
                                "start": None,
                                "end": None,
                                "source": "InterPro"
                            })

        except Exception as e:
            print(f"Domain error: {e}")

        return result

    def fetch_diseases(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Step 7: 获取疾病关联
        """
        result = []

        try:
            url = f"{self.uniprot_base}/{uniprot_id}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()

                if "comments" in data:
                    for comment in data["comments"]:
                        if comment["commentType"] == "DISEASE":
                            if "disease" in comment:
                                disease = comment["disease"]
                                name = disease.get("description", "")
                                references = disease.get("references", [])
                                evidence = "literature" if references else "database"

                                result.append({
                                    "name": name,
                                    "evidence": evidence,
                                    "source": "UniProt"
                                })

        except Exception as e:
            print(f"Disease error: {e}")

        return result

    def fetch_variants(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Step 8: 获取变异信息
        """
        result = []

        try:
            url = f"{self.uniprot_base}/{uniprot_id}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()

                if "features" in data:
                    for feature in data["features"]:
                        if feature["type"] == "SequenceVariant":
                            var_type = feature.get("variantType", "")
                            description = feature.get("description", "")

                            # 判断致病性
                            pathogenicity = None
                            if "consequence" in feature:
                                pathogenicity = feature["consequence"]

                            result.append({
                                "type": var_type,
                                "description": description,
                                "pathogenicity": pathogenicity,
                                "position": feature.get("location", {}).get("start", {}).get("value"),
                                "source": "UniProt"
                            })

        except Exception as e:
            print(f"Variants error: {e}")

        return result

    def fetch_structure(self, uniprot_id: str) -> Dict[str, Any]:
        """
        Step 9: 获取结构信息
        """
        result = {
            "has_experimental_structure": False,
            "has_alphafold_prediction": False,
            "pdb_ids": [],
            "alphafold_url": None,
            "source": None
        }

        try:
            url = f"{self.uniprot_base}/{uniprot_id}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()

                # 查找 PDB 交叉引用
                if "uniProtKBCrossReferences" in data:
                    for ref in data["uniProtKBCrossReferences"]:
                        if ref["database"] == "PDB":
                            result["pdb_ids"].append(ref["id"])
                            result["has_experimental_structure"] = True

                # AlphaFold 检测（AlphaFold DB 在 UniProt 也有引用）
                for ref in data.get("uniProtKBCrossReferences", []):
                    if ref["database"] == "AlphaFoldDB":
                        result["has_alphafold_prediction"] = True
                        result["alphafold_url"] = f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}"
                        break

        except Exception as e:
            print(f"Structure error: {e}")

        return result

    def annotate(self, input_id: str, species: str = "human") -> Dict[str, Any]:
        """
        完整注释流程
        """
        # Step 0: ID 规范化
        mapping = self.normalize_id(input_id, species)
        if not mapping["success"] or not mapping["uniprot_id"]:
            return {
                "success": False,
                "error": f"无法标准化 ID: {input_id}",
                "query": {
                    "input": input_id,
                    "species": species
                }
            }

        uniprot_id = mapping["uniprot_id"]

        # 依次获取各个模块
        result = {
            "success": True,
            "query": {
                "input": input_id,
                "input_type": mapping["input_type"],
                "species": species,
                "normalized_uniprot_id": uniprot_id,
                "normalized_gene_symbol": mapping["gene_symbol"],
                "normalized_ensembl_id": mapping["ensembl_id"]
            },
            "basic_info": self.fetch_basic_info(uniprot_id),
            "GO": self.fetch_go_annotations(uniprot_id),
            "pathways": self.fetch_pathways(uniprot_id),
            "interactions": self.fetch_interactions_string(uniprot_id, species),
            "domains": self.fetch_domains(uniprot_id),
            "disease": self.fetch_diseases(uniprot_id),
            "variants": self.fetch_variants(uniprot_id),
            "structure": self.fetch_structure(uniprot_id)
        }

        # 表达模块需要更多数据，先留空框架
        result["expression"] = {
            "source": "Not fetched",
            "tissues_high": [],
            "cell_types": []
        }

        return result


annotator = ProteinAnnotator()


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/annotate', methods=['POST'])
def annotate():
    data = request.get_json()
    input_id = data.get('input_id', '').strip()
    species = data.get('species', 'human')

    if not input_id:
        return jsonify({"success": False, "error": "请输入蛋白质 ID"})

    result = annotator.annotate(input_id, species)
    return jsonify(result)


@app.route('/api/health', methods=['GET'])
def health():
    return jsonify({"status": "ok", "service": "ProteinClaw"})


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
