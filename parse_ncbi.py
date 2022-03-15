import os, sys
import gzip

def main():
    data_dir = "../../data/raw"
    gene_info_file = os.path.join(data_dir, "Homo_sapiens.gene_info.gz")
    gene2pubmed_file = os.path.join(data_dir, "gene2pubmed.gz")
    gene2ensembl_file = os.path.join(data_dir, "gene2ensembl.gz")
    gene2go_file = os.path.join(data_dir, "gene2go.gz")

    geneid_to_name, name_to_geneid, geneid_to_synonyms, geneid_to_taxid = parse_ncbi.get_geneid_symbol_mapping(gene_info_file)
    geneid_to_pubmeds = parse_ncbi.get_geneid_to_pubmeds(gene2pubmed_file)
    geneid_to_ensembl, geneid_to_accession, _ = get_geneid_to_ensembl(gene2ensembl_file)
    geneid_to_go, geneid_to_go_to_evidences, geneid_to_go_to_pubmeds, _ = get_geneid_to_go(gene2go_file)

    return

def get_geneid_symbol_mapping(file_name):
    """
    Parses Homo_sapiens.gene_info from NCBI. 
    The file can be downloaded at: https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
    """
    geneid_to_name = {} # now contains only the official symbol
    name_to_geneid = {}
    geneid_to_synonyms = {}
    name_to_synonyms = {}
    geneid_to_taxid = {}

    f = gzip.open(file_name,'rb')
    first_line = f.readline()
    for line in f:
        words = line.decode().strip('\n').split('\t')
        if len(words) == 2:
            geneid, symbol = words
        else:
            taxid, geneid, symbol, locus, synonyms = words[:5]
            synonyms = synonyms.split("|")
        #print(taxid, geneid, symbol, locus, synonyms)

        # Insert TaxID. There can only be one taxID for GeneID. If not, error
        if geneid not in geneid_to_taxid:
            geneid_to_taxid[geneid] = taxid
        else:
            if geneid_to_taxid[geneid] != taxid:
                print('Different taxIDs for the GeneID: {}\nFirst taxID: {}  Second taxID: {}'.format(geneid, geneid_to_taxid[geneid], taxid))
                sys.exit(10)

        geneid = geneid.strip() # strip in case mal formatted input file
        symbol = symbol.strip()
        if geneid == "" or geneid == None or symbol == "" or symbol == None:
            continue
        
        #geneid_to_names.setdefault(geneid, set()).add(symbol) 
        geneid_to_name[geneid] = symbol
        for synonym in synonyms:
            geneid_to_synonyms.setdefault(geneid, set())
            geneid_to_synonyms[geneid].add(synonym)
            name_to_synonyms.setdefault(symbol)
            name_to_synonyms[symbol].add(synonym)
        
        for symbol in [symbol] + synonyms: # added for synonym parsing
            if symbol in name_to_geneid: 
                if int(geneid) >= int(name_to_geneid[symbol]):
                    continue
                #print("Multiple geneids", name_to_geneid[symbol], geneid, symbol)
            name_to_geneid[symbol] = geneid
    
    return geneid_to_name, name_to_geneid, geneid_to_synonyms, geneid_to_taxid


def get_geneid_to_pubmeds(file_name, tax_id = "9606"):
    """
    Parses gene2pubmed file from NCBI.
    The file can be downloaded at: https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz
    """
    geneid_to_pubmeds = {}
    f = gzip.open(file_name,'rb')
    first_line = f.readline()
    for line in f:
        tax, geneid, pubmed_id = line.decode().strip().split("\t")
        if tax != tax_id:
            continue
        geneid_to_pubmeds.setdefault(geneid, set()).add(pubmed_id)
    return geneid_to_pubmeds


def get_geneid_to_ensembl(file_name, id_type="geneid"):
    """
    Parses gene2ensembl file from NCBI.
    The file can be downloaded at: https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
    """
    geneid_to_ensembl = {}
    geneid_to_accession = {}
    geneid_to_taxid = {}

    gene2ensembl_fd = gzip.open(file_name,'rb')

    first_line = gene2ensembl_fd.readline().decode()

    # Obtain a dictionary: "field_name" => "position"
    fields_dict = obtain_header_fields(first_line[1:])
    #tax_id GeneID  Ensembl_gene_identifier RNA_nucleotide_accession.version    Ensembl_rna_identifier  protein_accession.version   Ensembl_protein_identifier

    for line in gene2ensembl_fd:

        # Split the line in fields
        fields = line.decode().strip().split("\t")

        # Obtain the fields of interest
        taxid = fields[ fields_dict['tax_id'] ]
        geneid = fields[ fields_dict['GeneID'] ]
        ensembl_gene = fields[ fields_dict['Ensembl_gene_identifier'] ]
        accession_rna = fields[ fields_dict['RNA_nucleotide_accession.version'] ]
        ensembl_rna = fields[ fields_dict['Ensembl_rna_identifier'] ]
        accession_prot = fields[ fields_dict['protein_accession.version'] ]
        ensembl_prot = fields[ fields_dict['Ensembl_protein_identifier'] ]

        # Insert TaxID. There can only be one taxID for GeneID. If not, error
        if geneid not in geneid_to_taxid:
            geneid_to_taxid[geneid] = taxid
        else:
            if geneid_to_taxid[geneid] != taxid:
                print('Different taxIDs for the GeneID: {}\nFirst taxID: {}  Second taxID: {}'.format(geneid, geneid_to_taxid[geneid], taxid))
                sys.exit(10)

        # Insert Ensembl
        geneid_to_ensembl.setdefault(geneid, {})

        if ensembl_gene != '-':
            ensembl_gene = ensembl_gene.split('.')[0]
            geneid_to_ensembl[geneid].setdefault('gene', [])
            geneid_to_ensembl[geneid]['gene'].append(ensembl_gene)
        if ensembl_rna != '-':
            ensembl_rna = ensembl_rna.split('.')[0]
            geneid_to_ensembl[geneid].setdefault('rna', [])
            geneid_to_ensembl[geneid]['rna'].append(ensembl_rna)
        if ensembl_prot != '-':
            ensembl_prot = ensembl_prot.split('.')[0]
            geneid_to_ensembl[geneid].setdefault('protein', [])
            geneid_to_ensembl[geneid]['protein'].append(ensembl_prot)

        # Insert Accession
        geneid_to_accession.setdefault(geneid, {})

        if accession_rna != '-':
            geneid_to_accession[geneid].setdefault('rna', [])
            geneid_to_accession[geneid]['rna'].append(accession_rna.split('.'))
        if accession_prot != '-':
            geneid_to_accession[geneid].setdefault('protein', [])
            geneid_to_accession[geneid]['protein'].append(accession_prot.split('.'))

    gene2ensembl_fd.close()

    return geneid_to_ensembl, geneid_to_accession, geneid_to_taxid


def get_geneid_to_go(file_name):
    """
    Parses gene2go file from NCBI.
    The file can be downloaded at: https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    """
    geneid_to_go = {}
    geneid_to_go_to_evidences = {}
    geneid_to_go_to_pubmeds = {}
    geneid_to_taxid = {}

    gene2go_fd = gzip.open(file_name,'rb')
    first_line = gene2go_fd.readline().decode()

    # Obtain a dictionary: "field_name" => "position"
    fields_dict = obtain_header_fields(first_line[1:])
    #tax_id GeneID  GO_ID   Evidence    Qualifier   GO_term PubMed  Category

    for line in gene2go_fd:

        # Split the line in fields
        fields = line.decode().strip().split("\t")

        # Obtain the fields of interest
        taxid = fields[ fields_dict['tax_id'] ]
        geneid = fields[ fields_dict['GeneID'] ]
        go_id = fields[ fields_dict['GO_ID'] ].upper()
        evidence = fields[ fields_dict['Evidence'] ].upper()
        qualifier = fields[ fields_dict['Qualifier'] ].lower()
        go_term = fields[ fields_dict['GO_term'] ].lower()
        pubmed = fields[ fields_dict['PubMed'] ]
        category = fields[ fields_dict['Category'] ].lower()

        # Insert TaxID. There can only be one taxID for GeneID. If not, error
        if geneid not in geneid_to_taxid:
            geneid_to_taxid[geneid] = taxid
        else:
            if geneid_to_taxid[geneid] != taxid:
                print('Different taxIDs for the GeneID: {}\nFirst taxID: {}  Second taxID: {}'.format(geneid, geneid_to_taxid[geneid], taxid))
                sys.exit(10)

        # Insert GO
        go_id = go_id.split('GO:')[1]
        geneid_to_go.setdefault(geneid, set())
        geneid_to_go[geneid].add(go_id)

        if evidence != '-':
            geneid_to_go_to_evidences.setdefault(geneid, {})
            geneid_to_go_to_evidences[geneid].setdefault(go_id, set())
            geneid_to_go_to_evidences[geneid][go_id].add(evidence)
        else:
            print('No evidence for GeneID {} and GO ID {}!'.format(geneid, go_id))
            sys.exit(10)

        if pubmed != '-':
            geneid_to_go_to_pubmeds.setdefault(geneid, {})
            geneid_to_go_to_pubmeds[geneid].setdefault(go_id, set())
            geneid_to_go_to_pubmeds[geneid][go_id].add(pubmed)

    gene2go_fd.close()

    return geneid_to_go, geneid_to_go_to_evidences, geneid_to_go_to_pubmeds, geneid_to_taxid



def get_unigene_to_geneids(file_name, prefix = "Hs."):
    """
    Parses gene2unigene file from NCBI.
    The file can be downloaded at: https://ftp.ncbi.nlm.nih.gov/gene/DATA/ARCHIVE/gene2unigene
    """
    unigene_to_geneids = {}
    with open(file_name) as f:
        f.readline()
        for line in f:
            geneid, unigene = line.strip().split("\t")
            if not unigene.startswith(prefix):
                continue
            unigene_to_geneids.setdefault(unigene, set()).add(geneid)
        #for unigene, geneids in unigene_to_geneids.iteritems():
        #        if len(geneids) > 1:
        #            print(unigene, geneids)
    return unigene_to_geneids


def get_homology_mapping(file_name, tax_id, from_tax_id="9606", symbol_type="geneid"):
    """
    Parses Homologene and obtains the homology mappings from one species to another.
    The file can be downloaded at: https://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data 
    File_name: Homologene data file
    tax_id: Tax id of the species to the genes of which the mapping 
    will be done (e.g., to mouse genes, from human genes)
    from_tax_id: Tax id from which the mapping will be done (default is human)
    symbol_type: geneid | symbol
    Tax ids for popular organisms
    Homo sapiens: 9606
    Mus musculus: 10090
    Rattus norvegicus: 10116
    Drosophila melanogaster: 7227
    Caenorhabditis elegans: 6239
    Saccharomyces cerevisiae: 4932
    Escherichia coli: 562
    Arabidopsis thaliana: 3702
    """
    # 3  9606  34  ACADM  4557231  NP_000007.1
    if symbol_type == "symbol":
        idx = 3
    elif symbol_type == "geneid":
        idx = 2
    else:
        raise ValueError("Unknown symbol type: {}".format(symbol_type))
    # Parse group info
    group_to_taxid_to_geneid = {}
    with open(file_name) as f:
        for line in f:
            words = line.strip("\n").split("\t")
            group, taxid = words[:2]
            geneid = words[idx]
            d = group_to_taxid_to_geneid.setdefault(group, {})
            d[taxid] = geneid
    # Get geneid mapping
    geneid_to_geneid = {}
    for group, taxid_to_geneid in group_to_taxid_to_geneid.iteritems():
        if tax_id not in taxid_to_geneid or from_tax_id not in taxid_to_geneid:
            continue
        geneid_to_geneid[taxid_to_geneid[from_tax_id]] = taxid_to_geneid[tax_id]
    return geneid_to_geneid, group_to_taxid_to_geneid 


def obtain_header_fields(first_line):
    """ 
    Obtain a dictionary: "field_name" => "position" 
    """
    fields_dict = {}

    header_fields = first_line.strip().split("\t")
    for x in range(0, len(header_fields)):
        fields_dict[header_fields[x]] = x

    return fields_dict


if __name__ == "__main__":
    main()