from requests import HTTPError
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from Bio import Entrez
# from langchain.llms import OpenAI
from langchain_community.llms import OpenAI
import time
import streamlit.components.v1 as components



## helper function to search pubmed
@st.cache_data
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

def stream_data(text):
    for word in text.split():
        yield word + " "
        time.sleep(0.02)

def generate_response(input_text):
    llm = OpenAI(temperature=0.7, openai_api_key=openai_api_key)
    output = llm(input_text)
    # st.info(llm(input_text))
    with st.chat_message("user"):
        # st.write("Hello Commander 👋")
        st.write_stream(stream_data(output))


# web app structure
st.title("Nutrigenomics One Stop Shop")

st.write("Welcome to the Nutrigenomics Dashboard! This is a demo of the dashboard for the HATCH2024 project by IOB UGA team.")

st.header("BiteSizedGenes")
st.write("To search for articles related to a gene of interest and summarize them using OpenAI's GPT-3.5 Turbo model.")

openai_api_key = st.text_input('Please input your own OpenAI API Key', type='password')

# Search box
keyword = st.text_input('Enter a keyword to search for articles:', '')
show_prompt = st.toggle('show prompt for OpenAI summarization', False)

if st.button("Use MTHFR as example input"):
    keyword = "MTHFR"

if keyword:
    with st.spinner(f"Searching for articles related to {keyword}..."):
    # st.write(f"Searching for articles related to {keyword}...")
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
            
            st.subheader("OpenAI summarization")

            # if st.button('Summarize'):
            if not openai_api_key.startswith('sk-'):
                st.warning('Please enter your own OpenAI API key!', icon='⚠')
            else:
                # overview = f"The DataFrame has {article_df.shape[0]} rows and {article_df.shape[1]} columns. The columns are: {', '.join(article_df.columns)}."
                overview = "you are an AI assistant on a spaceship helping human to understand the info of research on a dataset. the main mission is helping the human to learn about the gene and potential health related function or symptom the gene of interest would cause. make sure to greet the commander with some cheerful opening. The following stream of information are the latest 7 publication related to the keyword {keyword} from PubMed."
                content = article_df.to_string()
                prompt_text = f"{overview}\n\nhere are all the info in a dataframe: {content}\n\nBased on the above data, can you summarize them into a 150 words report for the astronaut commander in a nutrigenomics one stop shop?"

                generate_response(prompt_text)
                

                if show_prompt:
                    # st.success('Prompt generated successfully!')
                    st.info(f"Prompt used: \n\n{prompt_text}")
                
                    

        else:
            st.write("No articles found for the given keyword.")

st.header("Stellar SNPeek")
st.write("To check all the SNPs your astronauts have and compare with known pathogenic or disease-associated SNPs")
uploaded_file = st.file_uploader("Upload the VCF here")

# a button for example data
if st.button("Use example data"):
    uploaded_file = "./data/chr01/Astro1_MTHFR.vcf"

if uploaded_file is not None:

    st.write(f"you are viewing SNPs from {uploaded_file}")
    dataframe = pd.read_csv(uploaded_file, comment="#", sep="\t", header=None)
    dataframe.columns = ["#CHROM", "position", "ID", "reference", "alternative", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    dataframe = dataframe[["position", "reference", "alternative"]]
    st.write(dataframe)

    st.subheader("Known pathogenic or disease-associated SNPs.")
    Astr_SNPS_df = pd.read_csv("data/Astr_SNPS.csv")
    Astr_SNPS_df.columns = ["position", "Most Severe clinical significance", "condition"]
    st.write(Astr_SNPS_df)


st.header('FuelFolio Finder')

st.write("to find out the best food for your astronauts but also save the most fuel")

human_daily_nutrition_df = pd.read_csv("data/human_daily_nutrition.csv")
st.subheader("human daily nutrition")
st.dataframe(human_daily_nutrition_df)
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

st.write("for seed weight please search here: https://ser-sid.org/")

# embed streamlit docs in a streamlit app
components.iframe("https://ser-sid.org/", width=960, height=640, scrolling=True)

# footer
st.divider()
"Embark on a Journey of Taste and Health — Where Every Bite is an Adventure. © 2024. Dive deeper, pack smarter, and elevate your meals with Chromosome Culinary Crew. Bon Voyage on your culinary exploration!"

