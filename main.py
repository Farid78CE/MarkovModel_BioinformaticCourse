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
        # for debugging add /virtualEnvironment/src/
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
        nucleotidesProbabilityNew = {
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
        
        nucleotidesProbabilityNew['A'] = (nucleotides['A'] - 1) / (sequenceSize - 1) 
        nucleotidesProbabilityNew['T'] = (nucleotides['T'] - 1) / (sequenceSize - 1)
        nucleotidesProbabilityNew['C'] = (nucleotides['C'] - 1) / (sequenceSize - 1)
        nucleotidesProbabilityNew['G'] = (nucleotides['G'] - 1) / (sequenceSize - 1)

        return (nucleotides, nucleotidesProbability, nucleotidesProbabilityNew)

    def observedPairAlphabetFrequency(self, sequence:str):

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
        ObservedPairNucleotidesPossibility = self.calculateObservedProbabilityOfEachPattern(pairNucleotides, sequencePairsSize)
        return ObservedPairNucleotidesPossibility

    def countPattern(self, pairNucleotides: dict, sequence:str, patterns:list):
        for eachPattern in patterns:
            pairNucleotides[eachPattern] = sequence.count(eachPattern)
        return pairNucleotides

    def calculateObservedProbabilityOfEachPattern(self, pairNucleotides:dict, sequencePairsSize:int):
        for keys in pairNucleotides:
            pairNucleotides[keys] /= sequencePairsSize
        return pairNucleotides

    def calculateExpectedProbabilityOfEachPattern(self, nucleotideProbability:dict, observedPairAlphabetProbability:dict, pattern:list):
        # P(A|A) = p(AA)/p(A)
        expectedPairNucleotidesProbability = {
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

        for eachPattern in pattern:
            # P(C|A) = P(AC)/P(A)
            # 'AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT'
            denominatorPattern = eachPattern[0]
            expectedPairNucleotidesProbability[eachPattern] = observedPairAlphabetProbability[eachPattern] / nucleotideProbability[denominatorPattern]

        return expectedPairNucleotidesProbability

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
    # nucleotideProbabilityNew = tupleData[2]

    observedPairAlphabetProbability = mm.observedPairAlphabetFrequency(sequence)
    patterns = ['AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT']
    expectedPairAlphabetProbability = mm.calculateExpectedProbabilityOfEachPattern(nucleotideProbability, observedPairAlphabetProbability, patterns)
    fileStream = open(os.getcwd() + "/filestream.txt", 'a', encoding="utf8")
    fileStream.write(str(observedPairAlphabetProbability))
    fileStream.write(str(expectedPairAlphabetProbability))
    fileStream.close()
    for keys in observedPairAlphabetProbability:
        print(str(keys) + ":" + str(observedPairAlphabetProbability[keys]))
    print()
    for keys in expectedPairAlphabetProbability:
        print(str(keys) + ":" + str(expectedPairAlphabetProbability[keys]))