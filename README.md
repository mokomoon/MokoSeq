# MokoSeq

MokoSeq is a modular, workspace‑based bioinformatics platform built with Streamlit.  
It lets you upload sequence files, organize them into projects, assign categories, and preview them through a clean, interactive UI.

This project is still growing — it’s my personal space to learn, experiment, and build real bioinformatics tools as I improve.

---

## Features

- Workspace + project system  
- Upload FASTA / FASTQ files  
- Category assignment with automatic saving  
- Obsidian‑style sidebar file tree  
- Mobile‑friendly collapsible sidebar  
- Sequence preview (first 500 bp/aa)  
- Modular code structure for easy expansion  

---

## Getting Started

### Install dependencies

pip install -r requirements.txt

### Run the app

streamlit run app.py

### Your data folder


A `data/` folder will be created automatically — this is where your workspaces and projects live.

---

## Project Structure

MokoSeq/
│
├── app.py
├── requirements.txt
├── README.md
├── LICENSE
├── .gitignore
├── setup.py
├── pyproject.toml
│
├── mokoseq/
│   ├── analysis.py
│   ├── workspace.py
│   ├── categories.py
│   ├── filetree.py
│   ├── ui.py
│   └── utils.py
│
├── data/ (auto‑generated)
│
├── docs/
└── tests/


## Future Plans

---

## Why I Built This

I wanted a flexible, aesthetic, beginner‑friendly environment to explore bioinformatics.  
MokoSeq is my first full project — and I’ll keep improving it as I grow.

---

##  Future Plans

- ORF detection  
- QC metrics  
- Motif scanning  
- Metadata extraction  
- Visualizations and plots  
- Exportable reports  
- Better file management  

---

## License

MIT License.

