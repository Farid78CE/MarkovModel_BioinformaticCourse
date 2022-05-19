import os
import typing 
from Bio import SeqIO
import glob
import time 

class MarkovModel:
    def readFile(self,):
        curr:str = ""
        sequences:list = [] 
        curr = os.getcwd()
        os.chdir(str(curr) + "/fasta/")
        curr = os.getcwd()
        listOfFastaFiles:list = [] 

        for file in glob.glob("*.fasta"):
            listOfFastaFiles.append(file)
        
        for eachFastaFile in listOfFastaFiles:
            gene = curr + "/" + eachFastaFile
            records = SeqIO.parse(gene, 'fasta')
            for record in records:
                sequences.append(record.seq)
        


        return sequences[0]

    def alaphabetFrequency(self, sequence:str):
        nucleotides = {
            "A":0,
            "T":0,
            "C":0,
            "G":0
        }
        nucleotidesProbability = {
            "A":0.0,
            "T":0.0,
            "C":0.0,
            "G":0.0
        }

        sequenceSize = self.sequenceSize(sequence)
        

        for alphabets in sequence:
            if alphabets == "A":
                nucleotides['A'] += 1
            elif alphabets == "T":
                nucleotides['T'] += 1
            elif alphabets == "C":
                nucleotides['C'] += 1 
            elif alphabets == "G":
                nucleotides["G"] += 1 
            else: 
                print("Invalid Alphabet")

        nucleotidesProbability['A'] = nucleotides['A'] / sequenceSize
        nucleotidesProbability['T'] = nucleotides['T'] / sequenceSize
        nucleotidesProbability['C'] = nucleotides['C'] / sequenceSize
        nucleotidesProbability['G'] = nucleotides['G'] / sequenceSize
        

        return (nucleotides, nucleotidesProbability)

    def pairAlphabetFrequency(self, sequence:str):

        pairNucleotides = {
            'AA':0,
            'AC':0,
            'AG':0,
            'AT':0,
            'CA':0,
            'CC':0,
            'CG':0,
            'CT':0,
            'GA':0,
            'GC':0,
            'GG':0,
            'GT':0,
            'TA':0,
            'TC':0,
            'TG':0,
            'TT':0
        }

        pairNucleotidesPossibility = {
            'AA':0.0,
            'AC':0.0,
            'AG':0.0,
            'AT':0.0,
            'CA':0.0,
            'CC':0.0,
            'CG':0.0,
            'CT':0.0,
            'GA':0.0,
            'GC':0.0,
            'GG':0.0,
            'GT':0.0,
            'TA':0.0,
            'TC':0.0,
            'TG':0.0,
            'TT':0.0
        }

        sequenceSize = self.sequenceSize(sequence)
        sequencePairsSize = self.sequencePairSize(sequenceSize)
        patterns= ['AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT']
        pairNucleotides = self.countPattern(pairNucleotides, sequence, patterns)
        pairNucleotidesPossibility = self.calculateProbabilityOfEachPattern(pairNucleotides, sequencePairsSize)
        print(pairNucleotidesPossibility)

    def countPattern(self, pairNucleotides: dict, sequence:str, patterns:list):
        for eachPattern in patterns:
            pairNucleotides[eachPattern] = sequence.count(eachPattern)
        return pairNucleotides

    def calculateProbabilityOfEachPattern(self, pairNucleotides:dict, sequencePairsSize:int):
        for keys in pairNucleotides:
            pairNucleotides[keys] /= sequencePairsSize
        return pairNucleotides

    def sequenceSize(self, sequence:str):
        sequenceSize = len(sequence)        
        return sequenceSize
  
    def sequencePairSize(self, sequenceSize):
        if sequenceSize%2==0:
            pass
        else: 
            sequenceSize -= 1
        return sequenceSize
    
if __name__ == '__main__':
    mm = MarkovModel()
    sequence = mm.readFile()
    tupleData = mm.alaphabetFrequency(sequence)
    nucleotideFrequency = tupleData[0]
    nucleotideProbability = tupleData[1]
    mm.pairAlphabetFrequency(sequence)

    