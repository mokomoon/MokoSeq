"""
MokoSeq - Streamlit Workspace
--------------------------------------------

This is the main Streamlit application file. It handles:

1. Workspace selection
2. Project selection
3. File uploads
4. Category assignment
5. Category persistence
6. Obsidian-style sidebar tree
7. Mobile-friendly collapsible sidebar
8. Nested tabs for file analysis
9. Integration with the mokoseq/ Python package

"""

import streamlit as st
from pathlib import Path

# Import your internal modules
from mokoseq.workspace import ensure_workspace, ensure_project
from mokoseq.categories import load_categories, save_categories
from mokoseq.filetree import render_sidebar
from mokoseq.utils import read_file
from mokoseq.analysis import preview_sequence



# 1. PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="MokoSeq",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("MokoSeq Workspace")



# 2. ROOT DIRECTORY FOR USER DATA
# ============================================================

ROOT = Path("data")
ROOT.mkdir(exist_ok=True)



# 3. SIDEBAR COLLAPSE TOGGLE  (menu)
# ============================================================

st.markdown("""
<style>
.sidebar-toggle {
    position: fixed;
    top: 10px;
    left: 10px;
    font-size: 26px;
    cursor: pointer;
    z-index: 9999;
    padding: 4px 10px;
    border-radius: 6px;
    background: rgba(0,0,0,0.15);
    backdrop-filter: blur(6px);
}
@media (max-width: 900px) {
    .sidebar-toggle {
        top: 6px;
        left: 6px;
        font-size: 24px;
    }
}
</style>
<div class='sidebar-toggle' onclick='window.parent.postMessage({"toggle_sidebar": true}, "*")'>≡</div>
""", unsafe_allow_html=True)

st.markdown("""
<script>
window.addEventListener("message", (event) => {
    if (event.data.toggle_sidebar) {
        const sidebar = window.parent.document.querySelector("section[data-testid='stSidebar']");
        if (sidebar.style.display === "none") {
            sidebar.style.display = "block";
        } else {
            sidebar.style.display = "none";
        }
    }
});
</script>
""", unsafe_allow_html=True)



# 4. WORKSPACE SELECTION
# ============================================================

st.subheader("Workspace")

existing_workspaces = [p.name for p in ROOT.iterdir() if p.is_dir()]
workspace = st.selectbox("Select workspace", existing_workspaces)
new_workspace = st.text_input("Create new workspace")

if new_workspace:
    workspace = new_workspace
    ensure_workspace(workspace)

if not workspace:
    st.warning("Please select or create a workspace.")
    st.stop()



# 5. PROJECT SELECTION
# ============================================================

st.subheader("Project")

workspace_path = ROOT / workspace
existing_projects = [p.name for p in workspace_path.iterdir() if p.is_dir()]

project = st.selectbox("Select project", existing_projects)
new_project = st.text_input("Create new project")

if new_project:
    project = new_project
    ensure_project(workspace, project)

if not project:
    st.warning("Please select or create a project.")
    st.stop()



# 6. LOADING CATEGORIES
# ============================================================

project_path = workspace_path / project
categories = load_categories(project_path)



# 7. FILE UPLOAD 
# ============================================================

st.subheader("Upload Files")

uploaded_files = st.file_uploader(
    "Upload FASTA/FASTQ files",
    accept_multiple_files=True
)

if uploaded_files:
    files_dir = ensure_project(workspace, project)
    for uf in uploaded_files:
        with open(files_dir / uf.name, "wb") as f:
            f.write(uf.getbuffer())

files_dir = project_path / "files"
files = sorted([f.name for f in files_dir.iterdir()]) if files_dir.exists() else []



# 8. CATEGORY ASSIGNMENT
# ============================================================

st.subheader("Assign Categories")

for fname in files:
    current_cat = None
    for cat, flist in categories.items():
        if fname in flist:
            current_cat = cat

    col1, col2 = st.columns([3, 2])
    with col1:
        st.write(f"**{fname}**")

    with col2:
        chosen = st.selectbox(
            f"Category for {fname}",
            list(categories.keys()) + ["+ New Category"],
            index=list(categories.keys()).index(current_cat) if current_cat in categories else len(categories),
            key=f"cat_{fname}"
        )

        if chosen == "+ New Category":
            new_cat = st.text_input(f"New category name for {fname}", key=f"newcat_{fname}")
            if new_cat:
                categories.setdefault(new_cat, [])
                chosen = new_cat

        for cat in categories:
            if fname in categories[cat] and cat != chosen:
                categories[cat].remove(fname)

        categories.setdefault(chosen, [])
        if fname not in categories[chosen]:
            categories[chosen].append(fname)

save_categories(project_path, categories)



# 9. SIDEBAR TREE
# ============================================================

with st.sidebar:
    st.header(f"📁 {workspace}")
    render_sidebar(ROOT)



# 10. MAIN PANEL — CATEGORY SECTIONS + TABS
# ============================================================

st.header(f"Project: {project}")

for cat, flist in categories.items():
    with st.expander(cat, expanded=False):

        if not flist:
            st.info("No files in this category yet.")
            continue

        tabs = st.tabs(flist)

        for i, fname in enumerate(flist):
            with tabs[i]:
                st.subheader(fname)

                file_path = files_dir / fname
                seq_data = read_file(file_path)

                st.code(preview_sequence(seq_data), language="text")

                st.info("Analysis functions will go here.")
