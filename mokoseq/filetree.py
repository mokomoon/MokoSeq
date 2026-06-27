import streamlit as st
from pathlib import Path

def render_sidebar(root):
    for workspace in root.iterdir():
        if workspace.is_dir():
            st.markdown(f"## {workspace.name}/")
            for project in workspace.iterdir():
                if project.is_dir():
                    st.markdown(F"- {project.name}/")
