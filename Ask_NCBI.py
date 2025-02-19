from pprint import pprint as pp
from requests import Session
from secrets_api import NCBI_key
import zipfile
import os
import re
import shutil
import glob
import numpy as np
import pandas as pd
import time
from itertools import islice
import matplotlib.pyplot as plt
import seaborn as sns

class AskNCBI:
    '''This class is used to interact with the NCBI API. It has methods to retrieve gene information,'''
    
    '''Constructor of the AskNCBI class. It sets all the atributes of the object.'''
    def __init__(self, key, page_size=20, path_to_download='/Users/maksimsgolubovics/Python_VScode/NCBI_zip/'):
        
        
        """
        This is the constructor of the AskNCBI: class. It sets all the atributes of the object.
        key: str
            The API key that was given by NCBI for the user.
        page_size: int
            The number of records to return in each page of a paginated response.
        """
        
        
        self.apiurl = 'https://api.ncbi.nlm.nih.gov/datasets/v2'
        self.headers =  {'accept': 'application/json', 'api-key': key }     #to communicate the api key
        self.session = Session()                                            #Session is a class from requests module. It helps optimize requests.
        self.session.headers.update(self.headers)
        self.pathway = path_to_download   #pathway to the directory where the zip files will be downloaded
        self.page_size = page_size
    
    '''Here are all private methods used only inside the class.'''
    def __complete_url_path(self, path_pam, url):
        
        
        """
        Complete the url path by adding the path parameters.
        Parameters
        ----------
        path_pam : str or list of str
            The path parameters to add to the url.
        url : str
            The url to add the path parameters too.
        Returns
        -------
        str
            The url with the path parameters added.
        """
        
        
        if isinstance(path_pam, str):
           url += path_pam
        else:
           for id in path_pam:
            url += id
            if id == path_pam[-1]:
                break
            if len(path_pam)>1:
                url += '%2C'  
        return url
    
    def __iterate_reports(self, r, y):
        
        
        """
        Iterate through the reports in a paginated response and add the gene_id to the list.

        Parameters
        ----------
        r : requests.Response
            The response object containing the paginated response.
        y : list
            The list to append the gene_id to.
        """


        for x in range(len(r.json()['reports'])): #.json() is a method that returns a dictionary
            data = r.json()['reports'][x].get('gene', np.nan).get('gene_id', np.nan) # .get() is used to fill with nan values, so that KeyError wouldn't be raised in describe_genes()
            y.append(data)

    def __pretty_name(self, x):
        
        """
        This function takes a string as input and returns a cleaned up version of the string.
        The function removes all commas, parentheses, single quotes, and square brackets from the string, 
        and replaces all spaces with underscores. This function is used to generate a pretty name for the 
        gene_id, so that it can be used as a filename without any problems.
        Parameters
        ----------
        x : str
            The string to clean up.
        Returns
        -------
        str
            The cleaned up string.
        """


        return re.sub(r"[,\(\)'\[\]]", '', str(x)).replace(' ', '_') #used it only once... but maybe I will use it again)

    def __which_fasta(self, param):
        
        
        """
        This function takes a string parameter and returns the corresponding FASTA file extension.
        Parameters
        ----------
        param : str
            The string parameter to be used to get the FASTA file extension.
        Returns
        -------
        str
            The FASTA file extension.
        Raises
        ------
        ValueError
            If the parameter does not match with any of the keys in the FASTA dictionary.
        """


        try:
         FASTA = {'FASTA_RNA':'.fna'} #add more extension later on
         value = FASTA[param]
         return value
        except:
            raise ValueError(f'No FASTA extension found!')

    def __zip_unzip(self, r, name, f):
        
        
        """
        This function takes a response object, a name, and a file extension as input parameters. 
        It downloads the zip file from the given response, unzips it, and removes the zip file.
        It also checks if the downloaded file exists and returns a message accordingly.
        If the file does not exist, it removes the entire directory created during the unzip process.
        Parameters
        ----------
        r : response
            The response object from the NCBI API.
        name : str
            The name of the directory where the zip file is to be downloaded.
        f : str
            The file extension of the file to be downloaded.
        Returns
        -------
        str
            A message indicating whether the download was successful or not.
        """


        name = self.__pretty_name(name)
        pathway_zip = self.pathway+name+'.zip'
        pathway_ex = self.pathway+name

        with open(pathway_zip, 'wb') as file:
            for chunk in r.iter_content(chunk_size=1024):
                file.write(chunk)

        with zipfile.ZipFile(pathway_zip, 'r') as zip_ref:
            zip_ref.extractall(pathway_ex)
            os.remove(pathway_zip)
        #Perhaps it makes sense to create its own method for this
        if bool(glob.glob(f'{pathway_ex+'/ncbi_dataset/data/'}/*{str(f)}')) == True: #be careful, cause if .faa does't exist, it will delete all .fna aswell.
            return print('Download succeed!!!')
        else:
            try:
                shutil.rmtree(pathway_ex)
            except PermissionError as e:
                print(f'PermissionError: {e}')
            except Exception as e:
                print(f"Error: {e}")
            return print('Download failed!!!')

    def __chunk_iterator(self, iterable, page_size):
        iterator = iter(iterable)
        while chunk := list(islice(iterator, page_size)):
            yield chunk
    
    def __describe_genes(self, ids, 
                        columns):
        

        """
        This function takes a list of gene IDs and a taxon as input parameters, 
        and returns a pandas DataFrame containing the description of each gene.
        Parameters
        ----------
        ids : list of str
            A list of gene IDs to retrieve the description for.
        taxon : str, optional
            The taxon to filter the results by. The default is 'human'.
        columns : list of str, optional
            The columns to include in the pandas DataFrame. The default is
            ['symbol', 'gene_id', 'synonyms', 'taxname', 'type', 'chromosomes', 
             'transcript_count', 'description'].
        Returns
        -------
        pandas.DataFrame
            A pandas DataFrame containing the description of each gene.
        """


        if isinstance(columns, str):
            columns = [columns]
        data = self.full_gene_sum(ids)
        df = pd.DataFrame(columns=columns)

        for _ in range(len(data)):                 #two for-loops were used, because the ingoing data is a list of .json() dictionaries
            for x in range(len(data[_]['reports'])):
                gene_data = data[_]['reports'][x].get('gene')  
                df.loc[len(df)] = [gene_data.get(f'{x}', np.nan) for x in columns]
        return df

    '''Here are all public methods that can be used outside the class.'''
    def full_gene_sum(self, ids, page_token='', sample=float('inf')):
        
        
        """
        This method takes a list of gene IDs and a taxon as input parameters,
        and returns a list of dictionaries containing the full gene summary for each gene.
        The method paginates through the results using the 'page_token' and 'page_size' parameters,
        and returns the results in a list of dictionaries.
        If the 'sample' parameter is set to a finite number, the method will stop paginating
        once it has reached the specified number of results.
        Parameters
        ----------
        ids : list of str
            A list of gene IDs to retrieve the full gene summary for.
        taxon : str, optional
            The taxon to filter the results by. The default is 'human'.
        page_token : str, optional
            The page token to use for pagination. The default is an empty string.
        sample : int, optional
            The maximum number of results to return. The default is infinity.
        Returns
        -------
        list of dict
            A list of dictionaries containing the full gene summary for each gene.
        """


        url = self.apiurl + '/gene/id/'
        parameters = {'page_size': f'{self.page_size}',
                       'page_token': f'{page_token}',}
        url_1 = self.__complete_url_path(ids, url)
        r = self.session.get(url_1, params=parameters, stream=True)
        z = 1
        y = []
        y.append(r.json())

        while 'next_page_token' in r.json() and z<sample:
            parameters['page_token'] = r.json()['next_page_token']
            time.sleep(2)                                               #slowes everything down, so you won't be kicked out... needs some experimentation
            r = self.session.get(url_1, params=parameters, stream=True)
            page_token = r.json().get('next_page_token', np.nan)
            z+=1
            y.append(r.json())

        if len(y)>1:
             del y[-1]  # you need to delete the last entry because it is just count info and raises KeyError later
        return y
    
    def describe_genes(self,ids, columns=['symbol', 'gene_id', 'synonyms', 'taxname',
                                          'type', 'chromosomes', 'transcript_count', 
                                          'description']):
        df = pd.DataFrame()                                                           # Chuck size of 500 also needs experimantation
        for chunk in self.__chunk_iterator(ids, 500):                                 #it is necessary to chunk data input because it can be too much for URL or API (not sure) and at somepoint it kiks you out...
            df = pd.concat([df, self.__describe_genes(ids=chunk, columns=columns)], ignore_index=True) 
        return df

    def get_all_genes(self, taxon, page_token='', types='PROTEIN_CODING',sample=float('inf')):


        """
        This function takes a taxon and an optional page token as input parameters, 
        and returns a list of all gene IDs for that taxon.
        Parameters
        ----------
        taxon : str
            The taxon to filter the results by.
        page_token : str, optional
            The page token to use for pagination. The default is an empty string.
        types : str, optional
            The type of gene to filter the results by. The default is 'PROTEIN_CODING'.
        sample : int, optional
            The maximum number of results to return. The default is infinity.
        Returns
        -------
        list of str
            A list of all gene IDs for the given taxon.
        """


        url = self.apiurl + '/gene/taxon/' + f'{taxon}'
        z = 1
        y = []
        parameters = {'returned_content':'IDS_ONLY', 
                      'page_size':f'{self.page_size}', 
                      'page_token': f'{page_token}', 
                      'types': f'{types}'}
        r = self.session.get(url, params=parameters, stream=True)
        self.__iterate_reports(r,y)
        pp('Dowloading...')

        while 'next_page_token' in r.json() and z<sample:
            parameters['page_token'] = r.json()['next_page_token']
            r = self.session.get(url, params=parameters, stream=True)
            self.__iterate_reports(r,y)
            page_token = r.json().get('next_page_token', np.nan)
            z+=1
            time.sleep(0.5)                                         #slowes everything down, so you won't be kicked out... needs some experimentation

        return y

    def find_gene_ids(self, symbols, taxon='human'):


        """
        This function takes a list or a string of gene symbols and a taxon as input parameters, 
        and returns a list of the corresponding gene IDs for that taxon. But not in the order.
        Parameters
        ----------
        symbols : str/list of str
            The list of gene symbols to search for.
        taxon : str, optional
            The taxon to filter the results by. The default is 'human'.
        Returns
        -------
        list of str
            A list of the corresponding gene IDs for the given taxon.
        """


        url = self.apiurl + '/gene/symbol/'
        url_1 = self.__complete_url_path(symbols,url)
        url_1 += '/taxon/'+f'{taxon}'
        r = self.session.get(url_1)
        y = []
        self.__iterate_reports(r,y)
        return y
    
    def dowload_fasta(self, ids, taxon='human', param='FASTA_RNA'):


        """
        This function takes a list of gene IDs, a taxon and a parameter as input parameters, 
        and returns a list of the corresponding fasta sequences for that taxon.
        Parameters
        ----------
        ids : list of str
            The list of gene IDs to search for.
        taxon : str, optional
            The taxon to filter the results by. The default is 'human'.
        param : str, optional
            The type of fasta sequence to download. The default is 'FASTA_RNA'.
        Returns
        -------
        list of str
            A list of the corresponding fasta sequences for the given taxon.
        """


        url = self.apiurl + '/gene/id/'
        url_1 = self.__complete_url_path(ids, url)
        url_1 += '/download'
        parameters = {'include_annotation_type': param}
        r = self.session.get(url_1, params=parameters, stream=True)
        return self.__zip_unzip(r,ids, self.__which_fasta(param))
    


'''SHOWCASE'''
'''The exempale code to dowloand big amout of data'''
#ask_to = AskNCBI(key = NCBI_key, page_size=1000)
#human_genes = ask_to.get_all_genes(taxon='human',types='PROTEIN_CODING', sample=7)
#pp(human_genes)
#df = ask_to.describe_genes(ids=human_genes)
#pp(df)
#df.to_csv('path')

'''The example what can you do with data'''
#Prepare the dataset for work
#protein_genes = pd.read_csv('/Users/maksimsgolubovics/Python_VScode/csvData/human.csv')
#protein_genes['chromosomes'] = protein_genes['chromosomes'].str.strip('[]').str.strip("'").str.replace("'", '')#Cleaning the column for easier use
#protein_genes['chromosomes'] = protein_genes['chromosomes'].astype('category')
#category_order = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','X, Y','MT','Un'] #Creating an order for better visualization
#protein_genes['chromosomes'] = pd.Categorical(protein_genes['chromosomes'], categories=category_order, ordered=True)

#looking for interesting datapoints
#print(protein_genes[protein_genes['chromosomes'] == 'MT'])
#print(protein_genes[protein_genes['chromosomes'] == 'X, Y'])
#print(protein_genes[protein_genes['chromosomes'] == 'Un'])
#shows the amout of genes that by chromosomes and its transcript count
#grouped = protein_genes.groupby('chromosomes')['transcript_count'].value_counts()
#print(grouped)

#Visualization of gene distribution on chromosomes
#sns.set_palette('colorblind')
#sns.set_style("whitegrid", {"grid.color": "0.8", "grid.linestyle": ":"})
#sns.set_context('notebook')
#plt.xlabel('Chromosomes')
#plt.ylabel('Number of protein-coding genes')
#plt.title('The distribution of protein-coding genes on the chromosomes', fontweight="bold")
#sns.histplot(x='chromosomes',data=protein_genes)
#plt.show()

'''The example of downloading the fasta files'''
#ask_to = AskNCBI(key = NCBI_key, path_to_download='/Users/maksimsgolubovics/Python_VScode/NCBI_zip/')
#ask_to.dowload_fasta(ask_to.find_gene_ids(symbols='IL21'))
