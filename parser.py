""" 
Author: 
    Ashka Shah
Date Modified: 
    April 11, 2024
Description: 
    User defined python class specifically designed to parse through GTF files and phastCon files provided by 
    University of California Santa Cruz (UCSC)
"""
import pickle
import time
import pyBigWig

class phastCon:
    def __init__(self, filename):
        self.file = pyBigWig.open(filename)
        self.given = 0
        self.found = 0
        print("Please note that bigWig files have 0-based positioning.")
    
    def get_chromosomes(self):
        return self.file.chroms()
        
    def get_mean_conservation(self, chromosome: str, start: int, end: int):
        return self.file.stats(chromosome, start, end)[0]
    
    def get_mean_conservation_intervals(self, chromosome: str, intervals: list):
        scores = []
        for interval in intervals:
            if interval == (-1,-1):
                scores.append(-1)
            else:
                self.given += 1
                score = self.get_mean_conservation(chromosome, interval[0], interval[1])
                if score is None:
                    scores.append(-1)
                else:
                    self.found += 1
                    scores.append(score)
        return scores

class Transcript:
    def __init__(self, transcript_id, start, end, strand):
        self.transcript_id = transcript_id
        self.strand = strand
        self.start = start
        self.end = end     
        self.exons = []
        self.cds = []
        self.utr5 = []
        self.utr3 =[]
        self.introns = []
        self.exon_cons = []
        self.intron_cons = []
        self.cds_cons = []
        self.utr5_cons = []
        self.utr3_cons = []
        self.intronless = True
        self.avg_exon_cons = 0
        self.avg_cds_cons = 0
        self.avg_utr5_cons = 0
        self.avg_utr3_cons = 0
        
    def __len__(self):
        return abs(self.end - self.start)
    
    def __str__(self):
        return "Transcript Id: {0}\nStrand: {1}\nStart: {2}\nEnd: {3}\nExons: {4}\nIntrons: {5}\n5UTR: {6}\n3UTR: {7}".format(self.transcript_id, self.strand, 
                                                                                       self.start, self.end, self.exons, self.introns,self.utr5,self.utr3)
    def get_exon_lengths(self):
        lengths = []
        for exon in self.exons:
            lengths.append(exon[1]-exon[0])
        return lengths
    
    def get_cds_lengths(self):
        lengths = []
        for cds in self.cds:
            lengths.append(cds[1]-cds[0])
        return lengths
    
    def get_intron_lengths(self):
        lengths = []
        for intron in self.introns:
            lengths.append(intron[1]-intron[0])
        return lengths
    
    def add_exon(self, start: int, end: int, num: int):           
        if len(self.exons) + 1 == num:
            self.exons.append((start-1, end)) # converting to 0-based
        else:
            print("Exons missing from position {} to {} for: ".format(len(self.exons),num-1), self)
            for i in range(num - len(self.exons) - 1):
                self.exons.append((-1,-1)) # Information about the exon is not known
            self.exons.append((start-1, end))
            
    def add_CDS(self, start: int, end: int, num: int):
        if len(self.cds) + 1 == num:
            self.cds.append((start-1, end)) # converting to 0-based
        else:
            for i in range(num - len(self.cds) - 1):
                self.cds.append((-1,-1)) # Information about the exon is not known
            self.cds.append((start-1, end))
            
    def add_UTR5(self, start: int, end: int):
        self.utr5.append((start-1,end))
    
    def add_UTR3(self, start: int, end: int):
        self.utr3.append((start-1,end))
    
    def infer_introns(self):
        # Inference will be 0-based as well
        if len(self.exons) > 1:
            self.intronless = False
            for i in range(len(self.exons)-1):
                if self.exons[i] == (-1,-1) or self.exons[i+1] == (-1,-1):
                    self.introns.append((-1,-1)) # Info about the intron is not known because info about exon is missing
                else:
                    # First exon position is smaller than second exon position
                    if self.exons[i][0] < self.exons[i+1][0]:
                        if self.exons[i][1] == self.exons[i+1][0]:
                            self.introns.append((-1,-1))
                            print("Two adjacent exons found for transcript ID: {}".format(self.transcript_id))
                        else:
                            self.introns.append((self.exons[i][1], self.exons[i+1][0]))
                    # Second exon position is greater than first exon position
                    else:
                        if self.exons[i+1][1] == self.exons[i][0]:
                            self.introns.append((-1,-1))
                            print("Two adjacent exons found for transcript ID: {}".format(self.transcript_id))
                            
                        else:
                            self.introns.append((self.exons[i+1][1], self.exons[i][0]))
                        
    def set_intron_cons(self, cons: list):
        self.intron_cons = cons
        
    def set_cds_cons(self, cons: list):
        self.cds_cons = cons
        
    def set_utr5_cons(self, cons: list):
        self.utr5_cons = cons
        
    def set_utr3_cons(self, cons: list):
        self.utr3_cons = cons
    
    def set_avg_utr5_cons(self):
        utr_lens = []
        total = 0
        for utr in self.utr5:
            utr_lens.append(utr[1]-utr[0])
            total += utr[1]-utr[0]
        
        for i in range(len(utr_lens)):
            if self.utr5_cons[i] != -1:
                self.avg_utr5_cons += (utr_lens[i]/total) * self.utr5_cons[i]  
        
    def set_avg_utr3_cons(self):
        utr_lens = []
        total = 0
        for utr in self.utr3:
            utr_lens.append(utr[1]-utr[0])
            total += utr[1]-utr[0]
        
        for i in range(len(utr_lens)):
            if self.utr3_cons[i] != -1:
                self.avg_utr3_cons += (utr_lens[i]/total) * self.utr3_cons[i]    
    
    def set_avg_cds_cons(self):
        cds_lens = self.get_cds_lengths()
        total = 0
        for i in range(len(cds_lens)):
            if self.cds_cons[i] != -1:
                total += cds_lens[i]
        
        for i in range(len(cds_lens)):
            if self.cds_cons[i] != -1:
                self.avg_cds_cons += (cds_lens[i]/total) * self.cds_cons[i]         

class Gene:

    def __init__(self, gene_id: str, chromosome: str, gene_name: str):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chromosome = chromosome
        self.transcripts = {}
        self.is_splicable = False
        
    def __str__(self):
        return 'Gene Name: {0}\nGene Id: {1}\nChromosome: {2}\nSplice Variants: {3}\n'.format(self.gene_name, self.gene_id, 
                                                                                       self.chromosome, self.is_splicable)
    def __iter__(self):
        return iter(self.transcripts.values())
    
    def __len__(self):
        return len(self.transcripts)
    
    def transcript_exists(self, transcript_id: str):
        if transcript_id in self.transcripts:
            return True
        return False
    
    def get_transcript(self, transcript_id: str):
        if self.transcript_exists(transcript_id):
            return self.transcripts[transcript_id]
        else:
            return None
                 
    def add_transcript(self, transcript_id, start, end, strand): 
        if self.transcript_exists(transcript_id):
            print("Please investigate duplicated transcript id: {}".format(transcript_id))
        else:
            self.transcripts[transcript_id] = Transcript(transcript_id, start-1, end, strand)
            if not self.is_splicable and len(self.transcripts) > 1:
                self.is_splicable = True

class Parser:
    def __init__(self, species:str, transcript_filename: str, exon_filename: str, phastCon_filename: str):
        self.species = species
        self.transcript_filename = transcript_filename
        self.exon_filename = exon_filename
        self.phastCon_filename = phastCon_filename
        self.phastCon = phastCon(self.phastCon_filename)
        self.genes = {}

    def __len__(self):
        return len(self.genes)
    
    def __iter__(self):
        return iter(self.genes.values())
    
    def __getstate__(self):
        state = self.__dict__.copy()
        del state['phastCon']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def gene_exists(self, gene_id: str):
        if gene_id in self.genes:
            return True
        return False
    
    def get_gene(self, gene_id: str):
        if self.gene_exists(gene_id):
            return self.genes[gene_id]
        else:
            return None
        
    def add_gene(self, gene_id: str, gene_name: str, chromosome: str):
        # Create a new gene if it doesn't exist
        if not self.gene_exists(gene_id):
            self.genes[gene_id] = Gene(gene_id,chromosome,gene_name)
        return self.get_gene(gene_id)
    
    def parse_files(self):
        # Adding the gene and transcript information 
        tfptr = open(self.transcript_filename,'r')
        lines = tfptr.readlines()
        for line in lines:
            columns = line.strip("\n").split("\t")
            gene_id = columns[8].strip("\"") # Gene Id
            gene_name = columns[10].strip("\"") # Gene Name
            chromosome = columns[0] # Chromosome
            tr_start = int(columns[3]) # Transcript start
            tr_end = int(columns[4]) # Transcript end
            tr_strand = columns[6] # Transcript strand 
            tr_id = columns[9].strip("\"") # Transcript Id
            
            if tr_strand != ".":
                # Add the gene
                gene = self.add_gene(gene_id,gene_name,chromosome)  
                # Add the transcript to the gene
                gene.add_transcript(tr_id,tr_start,tr_end,tr_strand)
            else:
                print("Transcript with id-{} does not have strand information. So it will be excluded.")
                
        tfptr.close()
        
        # Adding the exon information 
        efptr = open(self.exon_filename,'r')
        lines = efptr.readlines()
        for line in lines:
            columns = line.strip("\n").split("\t")
            gene_id = columns[8].strip("\"") # Gene Id
            ex_start = int(columns[3]) # Transcript start
            ex_end = int(columns[4]) # Transcript end
            ex_num = int(columns[10]) # Exon number
            tr_strand = columns[6] # Transcript strand 
            tr_id = columns[9].strip("\"") # Transcript Id
            typ = columns[2].strip()
                
            if tr_strand != ".":
                gene = self.get_gene(gene_id)
                if gene:
                    transcript = gene.get_transcript(tr_id)
                    if transcript:
                        if typ == "exon":
                            transcript.add_exon(ex_start,ex_end,ex_num)
                        elif typ == "CDS":
                            transcript.add_CDS(ex_start,ex_end,ex_num)
                        elif typ == "5UTR":
                            transcript.add_UTR5(ex_start,ex_end)
                        else:
                            transcript.add_UTR3(ex_start,ex_end)
                    else:
                        print("Transcript with id-{} not found".format(tr_id))
                else:
                    print("Gene with id-{} not found".format(gene_id))
            else:
                print("Transcript with id-{} does not have strand information. So it will be excluded.")
                
        efptr.close()
        
    def populate_introns_conservation_scores(self):
        for gene in self:
            for transcript in gene:
                transcript.infer_introns()
                transcript.set_intron_cons(self.phastCon.get_mean_conservation_intervals(gene.chromosome,transcript.introns))
                transcript.set_cds_cons(self.phastCon.get_mean_conservation_intervals(gene.chromosome,transcript.cds))
                transcript.set_utr5_cons(self.phastCon.get_mean_conservation_intervals(gene.chromosome,transcript.utr5))
                transcript.set_utr3_cons(self.phastCon.get_mean_conservation_intervals(gene.chromosome,transcript.utr3))
                transcript.set_avg_cds_cons()
                transcript.set_avg_utr5_cons()
                transcript.set_avg_utr3_cons()
        print("Found {} out of {} requested conservation scores".format(self.phastCon.found,self.phastCon.given))
        print("Done parsing.")        