import streamlit as st
import matplotlib.pyplot as plt
import networkx as nx
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os

# Apply custom CSS for better styling
st.markdown("""
    <style>
        .main {
            background-color: #f0f8ff;
        }
        h1 {
            color: #004d99;
            font-size: 36px;
            font-weight: bold;
            text-align: center;
            text-shadow: 2px 2px 4px #99ccff;
        }
        h2, h3 {
            color: #003366;
            text-align: center;
        }
        .stTextArea, .stButton > button {
            border-radius: 10px;
            padding: 12px;
        }
        .stButton > button {
            background: linear-gradient(to right, #4CAF50, #2E8B57);
            color: white;
            font-weight: bold;
            border: none;
            transition: 0.3s;
            font-size: 18px;
            box-shadow: 3px 3px 6px rgba(0, 0, 0, 0.2);
        }
        .stButton > button:hover {
            background: linear-gradient(to right, #2E8B57, #1E5D3F);
            transform: scale(1.05);
        }
    </style>
""", unsafe_allow_html=True)

# Function to save user sequences to a FASTA file
def save_fasta(user_sequences):
    fasta_file = "sequences.fasta"
    with open(fasta_file, "w") as f:
        f.write(user_sequences)
    return fasta_file

# Function to build a phylogenetic tree
def build_phylogenetic_tree(file_path):
    try:
        alignment = AlignIO.read(file_path, "fasta")
        calculator = DistanceCalculator("identity")
        constructor = DistanceTreeConstructor(calculator, "upgma")
        phylo_tree = constructor.build_tree(alignment)

        tree_image_file = "tree.png"
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        Phylo.draw(phylo_tree, axes=ax)
        plt.savefig(tree_image_file)
        plt.close()

        return phylo_tree, tree_image_file
    except Exception as e:
        return None, str(e)

# Function to visualize SNP interaction network
def visualize_interaction_network(edges):
    G = nx.Graph()
    G.add_edges_from(edges)
    pos = nx.spring_layout(G)
    fig, ax = plt.subplots(figsize=(6, 4))
    nx.draw(G, pos, with_labels=True, node_color="#ff9900", edge_color="gray",
            node_size=2500, font_size=10, font_weight='bold', ax=ax)
    return fig

# Streamlit UI
st.title("🧬 Phylogenetic Tree & Sequence Analysis")

# User input for FASTA sequences
st.subheader("Enter your sequences in FASTA format:")
user_input = st.text_area("Example format:",
                          ">Seq1\nATGCGTACGTTAGTAACTG\n>Seq2\nATGCGTACGTTAGTACCTG\n>Seq3\nATGCGTACGTTGGTAACTG\n>Seq4\nATGCGGACGTTAGTAACTG",
                          height=200)

if st.button("🚀 Generate Phylogenetic Tree & Analyze Sequences"):
    if user_input.strip() == "":
        st.error("⚠️ Please enter valid sequences in FASTA format.")
    else:
        fasta_file = save_fasta(user_input)
        st.success("✅ Phylogenetic Tree Generated Successfully!")

        # Build tree
        tree, tree_image = build_phylogenetic_tree(fasta_file)
        if tree and os.path.exists(tree_image):
            st.image(tree_image, caption="🌳 Phylogenetic Tree", use_column_width=True)
        else:
            st.error(f"❌ Error generating tree: {tree_image}")

        # Visualize SNP interaction network
        st.subheader("🔗 SNP Interaction Network")
        example_edges = [("Seq1", "Seq2"), ("Seq2", "Seq3"), ("Seq3", "Seq4")]
        fig_network = visualize_interaction_network(example_edges)
        st.pyplot(fig_network)

# Footer
st.markdown("""
    <hr>
    <p style='text-align: center; font-size: 16px; color: #555;'>Thank You for visiting BIOPHYLO</p>
""", unsafe_allow_html=True)
