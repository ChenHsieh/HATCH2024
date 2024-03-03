from requests import HTTPError
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from Bio import Entrez

## helper function to search pubmed
@st.cache
def search_pubmed(keyword):
    # Your email here
    Entrez.email = "hatch2024@gmail.com"

    # Define the search term ("MTHFR" in this case) and the database
    search_term = keyword
    database = "pubmed"

    # Use Entrez.esearch to search for the keyword in the database
    search_handle = Entrez.esearch(db=database, term=search_term, retmax=10)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    # Extract the list of Ids for the found articles
    id_list = search_results["IdList"]
    
    # Fetch details of the articles using Entrez.efetch
    fetch_handle = Entrez.efetch(db=database, id=",".join(id_list), retmode="xml")
    articles = Entrez.read(fetch_handle)
    fetch_handle.close()

    # Process and print out some information about each article
    for article in articles['PubmedArticle']:
        # Extract basic information
        article_title = article['MedlineCitation']['Article']['ArticleTitle']
        authors_list = article['MedlineCitation']['Article']['AuthorList']
        authors = ', '.join([auth['LastName'] + " " + auth.get('Initials','') for auth in authors_list if 'LastName' in auth])
        pub_date = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
        year = pub_date.get('Year', 'Unknown year')

        print(f"Title: {article_title}")
        print(f"Authors: {authors}")
        print(f"Publication Year: {year}")
        print("-" * 100)
# Function to search PubMed
def search_pubmed(keyword):
    Entrez.email = "your_email@example.com"  # Use a valid email
    try:
        handle = Entrez.esearch(db="pubmed", term=keyword, retmax=7)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except HTTPError as error:
        st.error(f"An error occurred: {error}")
        return []

# Function to fetch details of a single PubMed article by its ID
def fetch_article_details(article_id):
    try:
        handle = Entrez.efetch(db="pubmed", id=article_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        return records
    except HTTPError as error:
        st.error(f"An error occurred: {error}")
        return None

def get_article_info(article, article_id):
    # Base URL for PubMed articles
    base_url = "https://pubmed.ncbi.nlm.nih.gov/"
    # Constructing the full URL for the article
    url = f"{base_url}{article_id}"
    
    title = article['MedlineCitation']['Article']['ArticleTitle']
    # Linking the title with the article URL
    title_link = {url}
    
    try:
        authors_list = article['MedlineCitation']['Article']['AuthorList']
        authors = ', '.join([author['LastName'] + " " + author.get('ForeName', '') for author in authors_list if 'LastName' in author])
    except KeyError:  # Handle articles without authors listed
        authors = "No authors listed"
    
    pub_date = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
    year = pub_date.get('Year', 'Unknown year')
    
    abstract = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', ['No abstract available'])[0]
    
    return title, title_link, authors, year, abstract

def stream_df(article_df):
    yield article_df
# web app structure
st.title("Nutrigenomics One Stop Shop")

st.write("Welcome to the Nutrigenomics Dashboard! This is a demo of the dashboard for the HATCH2024 project by IOB UGA team.")

st.header("Genes of interest")

# Search box
keyword = st.text_input('Enter a keyword to search for articles:', '')

if st.button("Use MTHFR as example input"):
    keyword = "MTHFR"

if keyword:
    ids = search_pubmed(keyword)  # Assuming search_pubmed function is defined elsewhere
    if ids:
        st.write("You can click on the article titles to view more details.")
        st.write(f"Showing the top {len(ids)} articles")
        
        articles_info = []
        for article_id in ids:
            details = fetch_article_details(article_id)
            if details:
                article = details["PubmedArticle"][0]
                title, title_link, authors, year, abstract = get_article_info(article, article_id)
                articles_info.append([title, title_link, authors, year, abstract])
        
        # Create a DataFrame
        article_df = pd.DataFrame(articles_info, columns=['Title', 'link', 'Authors', 'Year', 'Abstract'])
        article_df.set_index("Year", inplace=True)

        
        st.write(article_df)

        st.write(f"For more search result please visit: https://www.ncbi.nlm.nih.gov/search/all/?term={keyword}")

    else:
        st.write("No articles found for the given keyword.")

st.header("Individual Astronaut's SNPs")

uploaded_file = st.file_uploader("Upload the VCF here")

# a button for example data
if st.button("Use example data"):
    uploaded_file = "./data/Astro1_MTHFR.vcf"

if uploaded_file is not None:

    st.write(f"you are viewing SNPs from {uploaded_file}")
    dataframe = pd.read_csv(uploaded_file, comment="#", sep="\t", header=None)
    dataframe.columns = ["#CHROM", "position", "ID", "reference", "alternative", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    dataframe = dataframe[["position", "reference", "alternative"]]
    st.write(dataframe)

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