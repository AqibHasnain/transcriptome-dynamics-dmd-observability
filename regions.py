from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

# Copyright(C) 2009 Iddo Friedberg & Ian MC Fleming
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment
# https://biopython.org/wiki/Intergenic_regions
# heavily modded from the biopython wiki, but the structure comes from the author's above!
def get_interregions(genbank_path,intergene_length=1):
    ''' The regions between genes are termed intergenic regions. This function mines every one of
    those regions present in a genbank (.gbff) file iff the region has length >= intergene_length.
    The output can be saved as a fasta file with the intergenic regions and descriptions such as 
    location, strand, and unique ID'''

    seq_record = next(SeqIO.parse(open(genbank_path), 'genbank'))
    genome_size = len(seq_record)
    cds_list = []
    # Loop over the genome file, get the CDS features on each of the strands
    for feature in seq_record.features:
        if feature.type == "CDS":
            mystart = feature.location.start.position
            myend = feature.location.end.position
            strand = feature.strand
            if feature.strand == -1 or feature.strand == 1:
                cds_list.append((mystart,myend,strand))
            else:
                print("No strand indicated %d-%d. Assuming +\n" % (mystart, myend))
                cds_list.append((mystart, myend, 1))
    intergenic_records = []
    for i, pospair in enumerate(cds_list[1:]):
        # Compare current start position to previous end position
        # first take care of the case where the gene starts at position 0 
        if cds_list[i][0] == 0: # if start of cds is at position 0
            last_end = cds_list[-1][1]
            this_start = genome_size
            if this_start - last_end >= intergene_length:
                intergene_seq = seq_record.seq[last_end:this_start]
                intergenic_records.append(SeqRecord(intergene_seq,id="%s-ign-%d"%(seq_record.name,len(cds_list)),\
                        description="%d-%d"%(last_end+1,this_start),))
        # the rest of the regions 
        else:
            last_end = cds_list[i][1]
            this_start = pospair[0]
            if this_start - last_end >= intergene_length:
                intergene_seq = seq_record.seq[last_end:this_start]
                intergenic_records.append(SeqRecord(intergene_seq,id="%s-ign-%d"%(seq_record.name,i),\
                            description="%d-%d"%(last_end+1,this_start),))

    return intergenic_records

def get_promoterregions(datatable,intergenic_records,genbank_path,rel_dist_thresh=5000):
    '''datatable is a class object that contains gene names, locations, etc. instantiated with Datatable
    intergenic_records is essentially a fasta with the intergenic sequences (where promoters live) and 
    it is this functions job to, for every gene and its location in the datatable, grab the corresponding 
    promoter region from the intergenic_records. Can write the promoter sequences to a fasta after. Each 
    promoter sequences will have a unique id corresponding to the gene it was corresponding to. '''
    
    # first get the genome size
    genome_size = len(next(SeqIO.parse(open(genbank_path),'genbank')))
    
    promoter_records = []
    for i,location in enumerate(datatable.locations):
        
        strand = location[2] # '1' or '-1' for sense and antisense
        dist = np.inf
        
        if strand == 1: # only look upstream for intergenic region
            gene_start = location[0]
            for feature in intergenic_records:
                intergenic_location = [int(x) for x in feature.description.split('-')]
                intergenic_end = intergenic_location[1]
                new_dist = gene_start - intergenic_end
                rel_dist = genome_size - intergenic_end + gene_start
                if new_dist < 0: # Case 1 for strand 1, upstream constraint needs to be enforced using rel_dist
                    if rel_dist < rel_dist_thresh: # rel_dist b/w intergenic_end & gene_start within rel_dist_thresh to enforce upstream constraint  
                        if rel_dist < dist:
                            dist = rel_dist
                            promoter_seq = feature.seq
                            ign_start = intergenic_location[0]
                            ign_end = intergenic_location[1]
                elif new_dist >= 0: # Case 2 for strand 1, upstream constraint is naturally satisfied
                    if new_dist < dist:
                        dist = new_dist
                        promoter_seq = feature.seq
                        ign_start = intergenic_location[0]
                        ign_end = intergenic_location[1]   
            
        elif strand == -1: # only look downstream for intergenic region
            gene_end = location[1]
            for feature in intergenic_records:
                intergenic_location = [int(x) for x in feature.description.split('-')]
                intergenic_start = intergenic_location[0]
                new_dist = intergenic_start - gene_end
                rel_dist = genome_size - gene_end + intergenic_start
                if new_dist >= 0: # Case 1 for strand -1, downstream constraint is naturally satisfied
                    if new_dist < dist:
                        dist = new_dist
                        promoter_seq = feature.seq.reverse_complement()
                        ign_start = intergenic_location[1]
                        ign_end = intergenic_location[0]
                elif new_dist < 0: # Case 2 for strand -1, downstream constraint enforced using rel_dist
                    if rel_dist < rel_dist_thresh:
                        if rel_dist < dist:
                            dist = rel_dist
                            promoter_seq = feature.seq.reverse_complement()
                            ign_start = intergenic_location[1]
                            ign_end = intergenic_location[0]
        
        promoter_records.append(SeqRecord(promoter_seq,id="%s_%s"%(datatable.locus_tags[i],str(strand)),\
                    description="%d-%d"%(ign_start,ign_end),name="%s%s"%('Promoter sequence for ',datatable.locus_tags[i]),))

    return promoter_records


def get_promoterregions_v2(locations,locus_tags,intergenic_records,genbank_path,rel_dist_thresh=5000):
    '''datatable is a class object that contains gene names, locations, etc. instantiated with Datatable
    intergenic_records is essentially a fasta with the intergenic sequences (where promoters live) and 
    it is this functions job to, for every gene and its location in the datatable, grab the corresponding 
    promoter region from the intergenic_records. Can write the promoter sequences to a fasta after. Each 
    promoter sequences will have a unique id corresponding to the gene it was corresponding to. '''
    
    # first get the genome size
    genome_size = len(next(SeqIO.parse(open(genbank_path),'genbank')))
    
    promoter_records = []
    for i,location in enumerate(locations):
        
        strand = location[2] # '1' or '-1' for sense and antisense
        dist = np.inf
        
        if strand == 1: # only look upstream for intergenic region
            gene_start = location[0]
            for feature in intergenic_records:
                intergenic_location = [int(x) for x in feature.description.split('-')]
                intergenic_end = intergenic_location[1]
                new_dist = gene_start - intergenic_end
                rel_dist = genome_size - intergenic_end + gene_start
                if new_dist < 0: # Case 1 for strand 1, upstream constraint needs to be enforced using rel_dist
                    if rel_dist < rel_dist_thresh: # rel_dist b/w intergenic_end & gene_start within rel_dist_thresh to enforce upstream constraint  
                        if rel_dist < dist:
                            dist = rel_dist
                            promoter_seq = feature.seq
                            ign_start = intergenic_location[0]
                            ign_end = intergenic_location[1]
                elif new_dist >= 0: # Case 2 for strand 1, upstream constraint is naturally satisfied
                    if new_dist < dist:
                        dist = new_dist
                        promoter_seq = feature.seq
                        ign_start = intergenic_location[0]
                        ign_end = intergenic_location[1]   
            
        elif strand == -1: # only look downstream for intergenic region
            gene_end = location[1]
            for feature in intergenic_records:
                intergenic_location = [int(x) for x in feature.description.split('-')]
                intergenic_start = intergenic_location[0]
                new_dist = intergenic_start - gene_end
                rel_dist = genome_size - gene_end + intergenic_start
                if new_dist >= 0: # Case 1 for strand -1, downstream constraint is naturally satisfied
                    if new_dist < dist:
                        dist = new_dist
                        promoter_seq = feature.seq.reverse_complement()
                        ign_start = intergenic_location[1]
                        ign_end = intergenic_location[0]
                elif new_dist < 0: # Case 2 for strand -1, downstream constraint enforced using rel_dist
                    if rel_dist < rel_dist_thresh:
                        if rel_dist < dist:
                            dist = rel_dist
                            promoter_seq = feature.seq.reverse_complement()
                            ign_start = intergenic_location[1]
                            ign_end = intergenic_location[0]
        
        promoter_records.append(SeqRecord(promoter_seq,id="%s_%s"%(locus_tags[i],str(strand)),\
                    description="%d-%d"%(ign_start,ign_end),name="%s%s"%('Promoter sequence for ',locus_tags[i]),))

    return promoter_records




















    

