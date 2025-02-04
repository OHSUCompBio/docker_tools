import re

from edu.ohsu.compbio.annovar.annovar_parser import AnnovarVariantFunction

class VariantTranscript(AnnovarVariantFunction):
    '''
    '''
    def __init__(self, chromosome, position, ref, alt):
        '''
        Create new VariantTranscript using chromosome, position, ref, and alt.
        '''
        super().__init__(chromosome, position, ref, alt)

        self.protein_transcript = None
        
        # An hgvs.sequencevariant.SequenceVariant with g. notation 
        self.sequence_variant = None
        
        # Reference sequence of 10 bases before and after the variant 
        self.reference_context = None

    def __eq__(self, obj):
        if obj is None:
            return False
        elif super().__eq__(obj):
            return self.protein_transcript == obj.protein_transcript
        else:           
            return False
    
    def __str__(self):
        '''
        String representation of the VariantTranscript object
        '''
        return f'[AnnovarVariantFunction: genotype={self.chromosome}-{self.position}-{self.reference}-{self.alt}-transcript={self.refseq_transcript}-{self.protein_transcript}' + ('-splicing' if self.splicing else '') 

    def get_copy(self):
        '''
        Return a copy of this trancript 
        ''' 
        transcript = VariantTranscript(self.chromosome, self.position, self.reference, self.alt)
        transcript.variant_effect = self._noneIfEmpty(self.variant_effect)
        transcript.variant_type = self._noneIfEmpty(self.variant_type)
        transcript.hgvs_amino_acid_position = self._noneIfEmpty(self.hgvs_amino_acid_position)
        transcript.hgvs_base_position = self._noneIfEmpty(self.hgvs_base_position)
        transcript.exon = self._noneIfEmpty(self.exon)
        transcript.hgnc_gene = self._noneIfEmpty(self.hgnc_gene)
        transcript.sequence_variant = self._noneIfEmpty(self.sequence_variant)
        transcript.hgvs_c_dot = self._noneIfEmpty(self.hgvs_c_dot)
        transcript.hgvs_p_dot_one = self._noneIfEmpty(self.hgvs_p_dot_one)
        transcript.hgvs_p_dot_three = self._noneIfEmpty(self.hgvs_p_dot_three)
        transcript.splicing = self._noneIfEmpty(self.splicing)
        transcript.refseq_transcript = self._noneIfEmpty(self.refseq_transcript)
        transcript.protein_transcript = self._noneIfEmpty(self.protein_transcript)
        transcript.reference_context = self._noneIfEmpty(self.reference_context)

        return transcript

    def _noneIfEmpty(self, value: str):
        '''
        Return None if the string is an empty string.
        ''' 
        if value == '':
            return None
        return value
    
    def __lt__(self, other):
        '''
        The overloaded ``__lt__`` function is used to compare variant transcripts. 
        
        Variants that have different genotypes are compared by concatenating several fields together and comparing the strings.
        Variants that have the same genotype and transcript are compared by counting how many non-null fields they have (see _get_self_score()).
            - When scores are the same then transcript accession and version are compared.
        '''
        # Break down the transcript accessions          
        p = re.compile("(NM_|CCDS)(\d+)\.(\d+)")
        
        self_genotype = f'{self.chromosome}-{self.position}-{self.reference}-{self.alt}'
        self_result = p.search(self.refseq_transcript)
        self_transcript_type = self_result.group(1)
        self_accession = int(self_result.group(2))
        self_version = int(self_result.group(3))

        other_genotype = f'{other.chromosome}-{other.position}-{other.reference}-{other.alt}'
        other_result = p.search(other.refseq_transcript)
        other_transcript_type = other_result.group(1)
        other_accession = int(other_result.group(2))
        other_version = int(other_result.group(3))
        
        # Scoring is not meant to be used for comparing transcripts of different genotypes. If the variants are different or the 
        #transcripts are of different type (refseq and ccds) then comparison is presumably for the purpose of sorting rather ranking, 
        # so we compare the string values rather than scores. 
        if self_genotype != other_genotype or self_transcript_type != other_transcript_type:
            left = f"{self_genotype} {self.refseq_transcript} {self.hgvs_c_dot} {self.protein_transcript} {self.hgvs_p_dot_three}"
            right = f"{other_genotype} {other.refseq_transcript} {other.hgvs_c_dot} {other.protein_transcript} {other.hgvs_p_dot_three}"
            return left < right 

        # Return comparison based on score        
        self_score = self._get_self_score()
        other_score = other._get_self_score()
        
        if self_score != other_score:
            return self_score < other_score

        # If the accessions are the same but with different versions then choose the one with the latest version.
        if self_accession == other_accession:
            return self_version < other_version
                
        # If the accessions are different then the older accession is better because it is assumed to be the most curated, so
        # the lower accession number has higher rank. 
        return self_accession > other_accession
        
    def _get_self_score(self):
        '''
        Return a score based on how many the fields are filled in. This method can be improved in the following ways
        - Give c_dot a higher score when the value comes from HGVS/UTA rather than Annovar
        - Scoring some fields (e.g. g-dot, protein protein script, exon) may not be relevant for non coding variants. 
        '''
        score = 0
        
        # cDNA Fields
        if self.refseq_transcript:
            score += 1

        if self.hgvs_c_dot:
            score += 1
            
        if self.variant_type:
            score += 1
            
        if self.hgvs_base_position:
            score += 1

        if self.exon:
            score += 1
        
        if self.hgnc_gene:
            score += 1
                
        # Protein fields        
        if self.protein_transcript and self.hgvs_amino_acid_position and self.hgvs_p_dot_one and self.hgvs_p_dot_three:
            score += 1
            
        if self.variant_effect:            
            score += 1

        return score