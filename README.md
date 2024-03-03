# HATCH2024 - Nutrigenomics One Stop Shop by Team IOB UGA

## Team Members

We are a group of PhD students from the Institute of Bioinformatics at University of Georgia. Teams members are:

Chen Hsieh
Darrian Royce Talamentes
Mary-Frances Behnke

## Project Description

We are building a suite for commanders who is in charge of making sure that the crews are well-fed and healthy. The suite will be able to take in the crews' genetic information and provide a list of food that the soldiers should eat to maintain a healthy lifestyle and also optimize the payload.

## Project Goals

1. IDENTIFY the nutrigenomic profile from the given individuals.​

* Retrieve annotations for given genes from external databases.​

* SNP calling – identify variants and polymorphisms that result in unique nutritional needs (might require lit review, but also think about how we can do this bioinformatically)​

2. MAXIMIZE nutrition.​

* Here’s where we want to retrieve data from FooDB. I think it would be best to scrape this data directly from the web to reduce our app size, but we can also subset the full database first (e.g. only include plant-based foods) as another way to reduce the storage needed to by our app.​

* Keep in mind that the challenge wants it to be vegetable-based (we need to ask if that includes plant-based foods like tofu?), but there may be some flexibility here depending on which genes were provided in the fastas, which might require some quick literature review! (i.e., did the genes have to do with how the individual metabolizes macros like proteins and carbohydrates, or other nutrients like vitamins and minerals?)​

3. MINIMIZE payload weight.​

* Dehydration is a great way to reduce cargo volume and weight (FooDB includes water weight so we can do some calculations there)​

* However, Vitamins A and C are compromised/destroyed during dehydration process so they will need special consideration. I am not suggesting we devise a way to send hydrated food to space (I checked, this is not super feasible) – a better solution would be if the report generated a “flag” or “warning” for Vitamins A and C that recommends non-food supplements to fulfill those nutrient requirements.​

4. REPORT optimal vegetables/nutritional plans (and possibly alternates).​

* One way to quantitatively identify and report optimal vegetables for the space mission is to make a table based on nutritional value, which is directly populated from FooDB. From there we can develop some simple system that scores all vegetables in the database according to nutritional value per unit dry weight and scales the importance of a given nutrient based on the nutritional needs of the input group (see next slide).​

* Then, the nutritionally optimal fruits/vegetables are plotted on an x-y axis, where one axis is dehydrated weight, and the other axis is “nutritional value.”​

## Getting started

There are 3 main components to make use of the deliverables. The first is the SnakeMake pipeline we have for reproducible, scalable, and efficient data processing to do the variant calling. The Sencond is the cleaned database from FOODB for optimized recommendation. The third is the dashboard that we have built using Streamlit for the commander to research on genes and do data visualization. 

please replace the conda env as you see fit, we provided a sample conda env `data/envionments.yml` for you to use but it is not tested yet.

```bash
cd dashboard
conda activate ~/opt/anaconda3/envs/streamlit
streamlit run app.py
```
