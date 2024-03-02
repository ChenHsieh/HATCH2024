import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

st.title("Nutrigenomics Dashboard")

st.write("Welcome to the Nutrigenomics Dashboard! This is a demo of the dashboard for the HATCH2024 project by IOB UGA team.")

st.header("Genes with functional annotation")

df = pd.read_csv("data/Food.csv")

st.dataframe(df)  # Same as st.write(df)

st.header("Individual SNPs")

# Parameters
num_individuals = 10  # Number of individuals
num_snps = 5  # Number of SNPs

# Generate mock SNP IDs (e.g., SNP1, SNP2, ...)
snp_ids = ['SNP' + str(i+1) for i in range(num_snps)]

# Generate mock data for individuals
# For simplicity, we'll use only 'AA', 'AT', and 'TT' genotypes
genotypes = ['AA', 'AT', 'TT']
data = {snp: np.random.choice(genotypes, num_individuals) for snp in snp_ids}

# Generate individual IDs
individual_ids = ['Ind' + str(i+1) for i in range(num_individuals)]

# Create the DataFrame
df = pd.DataFrame(data, index=individual_ids)

st.write(df)

st.header('Interactive Scatterplot of Nutritional Rank vs Dehydrated Weight')

# mock nutrition and dehydration data
df = pd.DataFrame({
    "Nutritional Rank": np.random.randint(1, 100, 50), # Random ranks between 1 and 100
    "Dehydrated Weight": np.random.uniform(5, 20, 50) # Random weights between 5 and 20 grams
})
df


# Plot using Plotly
fig = px.scatter(df, x="Nutritional Rank", y="Dehydrated Weight", title="Nutritional Rank vs Dehydrated Weight", 
                 labels={"Nutritional Rank": "Nutritional Rank", "Dehydrated Weight": "Dehydrated Weight (grams)"},
                 template="plotly_white")

st.plotly_chart(fig, use_container_width=True)