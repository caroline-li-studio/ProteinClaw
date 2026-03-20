# ProteinClaw

ProteinClaw = 把"一个蛋白的所有已知注释"系统化取回并结构化输出

## 功能

输入一个蛋白质 ID（UniProt Accession / Gene symbol / Ensembl ID），自动从多个公共数据库获取并整合所有注释信息：

- ✅ ID 规范化与映射（UniProt ↔ Gene ↔ Ensembl）
- ✅ 基础信息（UniProt）
- ✅ GO 注释（Gene Ontology）
- ✅ Pathway 注释（KEGG + Reactome）
- ✅ 蛋白质相互作用（STRING）
- ✅ 结构域注释（Pfam + InterPro）
- ✅ 表达数据（GTEx + Human Protein Atlas）
- ✅ 疾病关联（OMIM + DisGeNET + GWAS Catalog）
- ✅ 变异信息（ClinVar + UniProt）
- ✅ 结构信息（AlphaFold + PDB）

## 设计原则

1. **每条信息必须带来源**
2. **每条信息最好带证据等级**

## 输出格式

统一 JSON Schema：

```json
{
  "query": {
    "input": "SMOC1",
    "input_type": "gene_symbol",
    "species": "human",
    "normalized_id": "Q9H3U1"
  },
  "basic_info": {},
  "GO": [],
  "pathways": [],
  "interactions": [],
  "domains": [],
  "expression": {},
  "disease": [],
  "variants": [],
  "structure": {}
}
```

## 快速开始

```bash
pip install -r requirements.txt
python app.py
```

访问 http://localhost:5000

## 架构

- 后端：Flask
- 前端：Bootstrap + 3Dmol.js（结构可视化可选）
- 数据来源：多个公共数据库 REST API
- 部署：支持 Vercel

## License

MIT
