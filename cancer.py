import os
import re
import time
import requests
import pandas as pd
from Bio import Entrez
from bs4 import BeautifulSoup

def search_geo_accessions(keyword, retmax=1000, retstart=0):
    Entrez.email = "swetaraibms@gmail.com"  
    search_results = []
    handle = Entrez.esearch(db="gds", term=keyword, retmax=retmax, retstart=retstart)
    record = Entrez.read(handle)
    handle.close()
    id_list = record.get("IdList", [])
    return id_list, int(record.get("Count", 0))

def fetch_geo_accession_details(id_list):
    Entrez.email = "swetaraibms@gmail.com"
    accession_numbers = set() 
    
    if id_list:
        handle = Entrez.efetch(db="gds", id=id_list, rettype="full", retmode="text")
        data = handle.read()  
        handle.close()
        accession_pattern = re.compile(r"GSE\d{6,7}")  # this can be changed according to need like if we need GEO accession numbers for 6 digits or 7 accordingly.
        found_accessions = accession_pattern.findall(data)
        accession_numbers.update(found_accessions)
    else:
        id_list is None

    return list(accession_numbers)

def fetch_all_geo_accessions(keyword, max_results=1000):
    all_accessions = set()
    retstart = 0
    while True:
        geo_ids, total_count = search_geo_accessions(keyword, retmax=max_results, retstart=retstart)
        if geo_ids is None:
            break
        else:
            accessions = fetch_geo_accession_details(geo_ids)
            all_accessions.update(accessions)
            retstart += max_results
            if retstart >= total_count:
                break
    return list(all_accessions)

def fetch_all_geo_accessions(keyword, max_results=1000):
    all_accessions = set()
    retstart = 0
    while True:
        geo_ids, total_count = search_geo_accessions(keyword, retmax=max_results, retstart=retstart)
        if not geo_ids:
            break
        all_accessions.update(fetch_geo_accession_details(geo_ids))
        retstart += max_results
        if retstart >= total_count:
            break
    return list(all_accessions)

def save_accessions_to_sheet(accessions, directory, filename="geo_accessions.xlsx"):
    file_path = os.path.join(directory, filename)
    df = pd.DataFrame(accessions, columns=["GEO Accession Number"])
    df.to_excel(file_path, index=False)
    print(f"Data saved to {file_path}")

def fetch_pubmed_ids_from_geo(accession_list):
    base_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    metadata = []

    for accession in accession_list:
        url = f"{base_url}?acc={accession}&view=full"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')

        pubmed_ids = []
        for a in soup.find_all('a', href=True):
            if 'pubmed' in a['href']:
                pubmed_id = a['href'].split('=')[-1]
                pubmed_ids.append(pubmed_id)

        metadata.append({
            'accession': accession,
            'pubmed_ids': ', '.join(pubmed_ids)  
        })

    return metadata

def clean_pubmed_ids(file_path, output_path):
    df = pd.read_excel(file_path)
    df['Pubmed_ID'] = df['pubmed_ids'].str.extract(r'/pubmed/(\d+)$')
    df.to_excel(output_path, index=False)
    print(f"Cleaned PubMed IDs saved to {output_path}")
    return df
def fetch_pubmed_html(pubmed_id):
    url = f'https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/'
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    try:
        response = requests.get(url, headers=headers, timeout=10)
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching PubMed ID {pubmed_id}: {e}")
        return None


def search_nct_in_abstract(html):
    soup = BeautifulSoup(html, 'html.parser')

    abstract_div = soup.find('div', class_='abstract-content selected')
    if abstract_div:
        abstract_text = abstract_div.get_text()
        nct_match = re.search(r'NCT\d+', abstract_text)
        if nct_match:
            return nct_match.group(0)
        trial_reg_div = soup.find('div', class_='trial-registration')
        if trial_reg_div:
            trial_reg_text = trial_reg_div.get_text()
            nct_match = re.search(r'NCT\d+', trial_reg_text)
            if nct_match:
                result = nct_match.group(0)
            else:
                result = None
            return result
    else:
        return None

def process_pubmed_ids(input_file, output_file):
    df = pd.read_excel(input_file).dropna(subset=['Pubmed_ID'])
    df['Pubmed_ID'] = pd.to_numeric(df['Pubmed_ID'], errors='coerce').dropna().astype(int)
    pubmed_ids = df['Pubmed_ID']
    results = []
    for i, pubmed_id in enumerate(pubmed_ids, start=1):
        print(f"Processing {i}/{len(pubmed_ids)}: PubMed ID {pubmed_id}")
        a = fetch_pubmed_html(pubmed_id)
        if a:
            nct_number = search_nct_in_abstract(a)
            if nct_number:
                results.append({'Pubmed_ID': pubmed_id, 'NCT Number': nct_number})
                print(f"  Found NCT Number: {nct_number}")
            else:
                results.append({'Pubmed_ID': pubmed_id, 'NCT Number': 'NCT Not Found'})
                print(f"  No NCT Number found.")
        else:
            results.append({'Pubmed_ID': pubmed_id, 'NCT Number': 'NCT Not Found'})
            print(f"  Failed to fetch data for PubMed ID {pubmed_id}.")
        time.sleep(1)
    res1 = pd.DataFrame(results)
    res1.to_excel(output_file, index=False)
    print(f"Processed NCT numbers saved to {output_file}")

def filter_nct(input_file, output_file):
    data = pd.read_excel(input_file)
    match_data = data[data['NCT Number'].str.match(r'^NCT\d+$', na=False)]
    match_data.to_excel(output_file, index=False)
    print(f"Filtered data with valid NCT numbers saved to {output_file}")

dir = '/Users/swetarai/Downloads'
keyword = "Cancer" 
geo_accessions = fetch_all_geo_accessions(keyword)
save_accessions_to_sheet(geo_accessions, dir, filename="geo_accessions_cancer.xlsx")

input_geo = os.path.join(dir, "geo_accessions_cancer.xlsx")
geo_df = pd.read_excel(input_geo)
geo_accessions = geo_df['GEO Accession Number'].tolist()
geo_metadata = fetch_pubmed_ids_from_geo(geo_accessions)
metadata_file = os.path.join(dir, "geo_pubmed_metadata_oncogene_cancer.xlsx")
pd.DataFrame(geo_metadata).to_excel(metadata_file, index=False)

cleaned_file_data = os.path.join(dir, "cleaned_geo_pubmed_metadata_cancer.xlsx")
clean_pubmed_ids(metadata_file, cleaned_file_data)
nct_output_file = os.path.join(dir, "output_excel_file_cancer.xlsx")
process_pubmed_ids(cleaned_file_data, nct_output_file)
final_output_file = os.path.join(dir, "output_with_valid_nct_cancer.xlsx")
filter_nct(nct_output_file, final_output_file)