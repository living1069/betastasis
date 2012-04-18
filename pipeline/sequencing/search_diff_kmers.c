
#include <stdio.h>

int main(int argc, char** argv) {
	
	// FIXME: Make the code non-strand specific.
	
	unsigned int MAX_KMERS = 65536;
	unsigned int kmers_done = 0;
	char* kmers[MAX_KMERS];
	uint64_t kmer_reads[MAX_KMERS];
	char line[1024];
	char nucs[] = "ACGT";
	
	// Initialize with all 8-mers.
	for (int k = 0; k < 65536; k++) {
		char* kmer = malloc(9);
		int nuc = k;
		for (int n = 0; n < 8; n++) {
			kmer[n] = (nuc % 4);
			nuc /= 4;
		}
		kmer[8] = '\0';
		kmers[k] = kmer;
	}
	
	while (1) {
		for (int k = 1; k < argc; k++) {
			FILE* file = fopen(argv[k], "r");
			
			while (1) {
				if (fgets(line, 1024, file) == NULL)
					break;
				
				if (line[0] == '>')
					continue;
				
				int line_len = strlen(line);
				
				for (int c = 0; c < line_len; c++) {
					
				}
				
				strstr(line, 
			}
			
			fclose(file);
		}
	}
	
	return 0;
}


