import streamlit as st

def collapsible_section(title, files):
    with st.expander(title):
        if not files:
            st.info("No files in this category.")
            return
        tabs = st.tabs(files)
        return tabs
