''' packages '''
import urllib.parse
import urllib.request




''' 1 fasta a uniprotID'''




''' 2 Uniprot a PDB '''

# P06401, Q9P7Q4, Q9VVG4, P38401, P11433, Q9Y223, P65206

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'PDB_ID',
'format': 'tab',
'query': 'P40925 P40926 O43175 Q9UM73 P97793'
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
print(response.decode('utf-8'))



''' 3 uniprot a alphafold '''


''' 4 alphafold aconseguir B-factors '''


''' 5 B factors a partir de pdb ''''
