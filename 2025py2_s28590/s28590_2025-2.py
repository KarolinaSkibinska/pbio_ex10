from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt


class NCBIRetriever:
    def __init__(self,email,api_key):
        self.email=email
        self.api_key=api_key
        Entrez.email=email
        Entrez.api_key=api_key
        Entrez.tool='BioScriptEx10'

    def search_taxid(self,taxid,min_length=0,max_length=float('inf')):
        print(f"Searching for records with taxID:{taxid}")
        try:
            handle=Entrez.efetch(db="taxonomy",id=taxid,retmode="xml")
            records=Entrez.read(handle)
            organism_name=records[0]["ScientificName"]
            print(f"Organism:{organism_name}(TaxID:{taxid})")
            search_term=f"txid{taxid}[Organism]"
            handle=Entrez.esearch(db="nucleotide",term=search_term,usehistory="y")
            search_results=Entrez.read(handle)
            count=int(search_results["Count"])

            if count == 0:
                print(f"No records found for{organism_name}")
                return None

            print(f"Found{count}records")

            self.webenv=search_results["WebEnv"]
            self.query_key=search_results["QueryKey"]
            self.count=count
            self.records_to_fetch=[]
            self.min_length=min_length
            self.max_length=max_length

            return count

        except Exception as e:
            print(f"Error searching TaxID{taxid}:{e}")
            return None

    def fetch_records(self,start=0,max_records=10):
        if not hasattr(self,'webenv') or not hasattr(self,'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            batch_size=min(max_records, 500)
            handle=Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )
            records=SeqIO.parse(handle,"genbank")
            filtered_records=[]

            for record in records:
                seq_length=len(record.seq)
                if self.min_length<=seq_length<=self.max_length:
                    filtered_records.append((record.id,seq_length))

            return filtered_records

        except Exception as e:
            print(f"Error fetching records:{e}")
            return []


def generate_csv(records,output_filename):
    df=pd.DataFrame(records,columns=["Accession Number","Sequence Length"])
    df.to_csv(output_filename,index=False)
    print(f"CSV report generated:{output_filename}")


def generate_plot(records,plot_filename):
    records=sorted(records,key=lambda x: x[1],reverse=True)
    accession_numbers=[record[0] for record in records]
    lengths=[record[1] for record in records]
    plt.figure(figsize=(10, 6))
    plt.plot(accession_numbers,lengths,marker='o',linestyle='-',color='b')
    plt.xlabel('Accession Number')
    plt.ylabel('Sequence Length')
    plt.title('Sequence Lengths by Accession Number')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(plot_filename)
    print(f"Plot saved as:{plot_filename}")


def main():
    email=input("Enter your email address for NCBI: ")
    api_key=input("Enter your NCBI API key: ")
    retriever=NCBIRetriever(email, api_key)
    taxid=input("Enter taxonomic ID (taxid) of the organism: ")
    min_length=int(input("Enter minimum sequence length: "))
    max_length=int(input("Enter maximum sequence length: "))
    count=retriever.search_taxid(taxid, min_length, max_length)
    if not count:
        print("No records found. Exiting.")
        return
    print("\nFetching filtered records...")
    records=retriever.fetch_records(start=0, max_records=20)
    if not records:
        print("No filtered records found. Exiting.")
        return

    output_csv=f"taxid_{taxid}_filtered_records.csv"
    generate_csv(records,output_csv)
    output_plot=f"taxid_{taxid}_sequence_lengths.png"
    generate_plot(records, output_plot)


if __name__ == "__main__":
    main()
