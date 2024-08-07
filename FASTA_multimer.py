from fasta import FASTQ

MC_sequences = {"seq1" : "TKCGLMEYHGKKKTKCPQWC",
                "seq2" : "TKCECMEYHGKKKTKCKVWC",
                "seq3" : "TKYECMEYHGKKKTKDRVWC",
                "seq4" : "TKYECMIYYGKKKTKEEQWC",
                "seq5" : "TKYECPLYPDKKKTKREIWY",
                "seq6" : "TKYDIDVMDYKKKYYEEFEE",
                "seq7" : "TKYDIDCMDSKKKYYEEFEE",
                "seq8" : "TKYDIDNMDSKKKYYEEFEE",
                "seq9" : "TKYYIDAMGSKKKEYEEFEE",
                "seq10" : "SKQFCRDEVYKKRQCTKYCC",
            }

#checking the outputs
print(MC_sequences.items())

output_path = 'APO_006_fasta.txt'
output_file = open(output_path, 'w')

#APO A-I sequence
target_sequence = "MKAAVLTLAVLFLTGSQARHFWQQDEPPQSPWDRVKDLATVYVDVLKDSGRDYVSQFEGSALGKQLNLKLLDNWDSVTSTFSKLREQLGPVTQEFWDNLEKETEGLRQEMSKDLEEVKAKVQPYLDDFQKKWQEEMELYRQKVEPLRAELQEGARQKLHELQEKLSPLGEEMRDRARAHVDALRTHLAPYSDELRQRLAARLEALKENGGARLAEYHAKATEHLSTLSEKAKPALEDLRQGLLPVLESFKVSFLSALEEYTKKLNTQ"

#Iterating through the dictionary
for seq_id, seq in MC_sequences.items():
    identifier_line = ">" + seq_id + "\n"
    output_file.write(identifier_line)
    sequence_line = target_sequence + ":" + seq + "\n"
    output_file.write(sequence_line)

#Closing the file
output_file.close()