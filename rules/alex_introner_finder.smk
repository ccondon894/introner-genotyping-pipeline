from Bio import SeqIO
import os
import subprocess
import collections
from collections import defaultdict
import io

# Reference paths
deplete = "data/deplete-reference/20220622_1749_Jacobs.fasta"
intronerized = "data/intronerized-reference/MpusillaCCMP1545_228_v3.0.fa"
names = ["CCMP490", "RCC114", "RCC373", "RCC465", "RCC629", "RCC647", "RCC692", "RCC693", "RCC833", "RCC835", "RCC1614", "RCC1698", "RCC2482", "RCC3052"]

# Known introns paths
known_introns_intronerized = "data/CCMP1545.introns_from_rnaseq.fa"
known_introns_deplete = "data/RCC1749.introns_from_rnaseq.fa"
indels_from_pairwise_graph = 'data/mpusilla_perm_pg.vcf'

rule all:
    input:
        expand("output/6-strain-haplotypes/{strain}.{reference}.fasta", strain=names, reference=["deplete", "intronerized"]),
        "output/filtered.introners.from_intron.deplete.fa",
        "output/filtered.introners.from_intron.intronerized.fa",
        "output/filtered.introners.from_indels.deplete.fa",
        "output/filtered.introners.from_indels.intronerized.fa",
        "output/filtered.introners.from_denovo.deplete.fa",
        "output/filtered.introners.from_denovo.intronerized.fa",
        "output/filtered.introners.from_megagraph_indels.deplete.fa",
        "output/filtered.introners.from_megagraph_indels.intronerized.fa",
        "output/calls/filtered.introners.from_intron.deplete.fa",
        "output/calls/filtered.introners.from_intron.intronerized.fa",
        "output/calls/filtered.introners.from_indels.deplete.fa",
        "output/calls/filtered.introners.from_indels.intronerized.fa",
        "output/calls/filtered.introners.from_denovo.deplete.fa",
        "output/calls/filtered.introners.from_denovo.intronerized.fa",
        "output/calls/filtered.introners.from_megagraph_indels.deplete.fa",
        "output/calls/filtered.introners.from_megagraph_indels.intronerized.fa"

rule prepare_indels:
    input:
        indels = indels_from_pairwise_graph
    output:
        depl="data/insertions-in-deplete.fa",
        intr="data/insertions-in-intronerized.fa"
    run:
        with open(input.indels, 'r') as f, \
            open(output.depl, 'w') as depl, open(output.intr, 'w') as intr:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                scaf, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
                if len(ref) > len(alt) and len(ref) > 40:
                    intr.write(f">Ref.Insertion.Intronerized.{scaf}:{pos}\n{ref}\n")
                elif len(alt) > len(ref) and len(alt) > 40:
                    depl.write(f">Ref.Insertion.Deplete.{scaf}:{pos}\n{alt}\n")

rule local_denovo:
    input:
        ref_fasta=lambda wildcards: deplete if wildcards.reference == "deplete" else intronerized
    output:
        vcf="output/1-find-nonref-insertions/{strain}/{reference}/assembled_ins.vcf.gz"
    threads: 60
    shell:
        """
        mkdir -p output/1-find-nonref-insertions/{wildcards.strain}/{wildcards.reference}
        python3 ext/INSurVeyor/insurveyor.py --threads {threads} --min-insertion-size 40 --min-stable-mapq 50 data/strains-bams/{wildcards.strain}.{wildcards.reference}.bam output/1-find-nonref-insertions/{wildcards.strain}/{wildcards.reference} {input.ref_fasta}
        """

rule vcf2fa:
    input:
        vcf="output/1-find-nonref-insertions/{strain}/{reference}/assembled_ins.vcf.gz"
    output:
        fasta="output/1-find-nonref-insertions/{strain}/{reference}/denovo.fa"
    shell:
        """
        python3 src/vcf2fa.py {input.vcf} {output.fasta}
        """
rule cat_fa:
    input:
        fa= lambda wildcards: expand("output/1-find-nonref-insertions/{strain}/intronerized/denovo.fa", strain=names) if wildcards.reference == "intronerized" else expand("output/1-find-nonref-insertions/{strain}/deplete/denovo.fa", strain=names)
    output:
        all_fa="output/1-find-nonref-insertions/{reference}.all.fa"
    shell:
        """
        cat {input.fa} > {output.all_fa}
        """

rule cluster_sequences:
    input:
        ref_insert_intr="data/insertions-in-intronerized.fa",
        ref_insert_depl="data/insertions-in-deplete.fa",
        ref_intron_intr=known_introns_intronerized,
        ref_intron_depl=known_introns_deplete,
        new_ins_intr=expand("output/1-find-nonref-insertions/{strain}/intronerized/denovo.fa", strain=names),
        new_ins_depl=expand("output/1-find-nonref-insertions/{strain}/deplete/denovo.fa", strain=names)
    output:
        fa="output/tmp/all.filt.fa",
        clusters="output/2-cluster-sequences/all.fa.clstr"
    threads: 120
    shell:
        """
        cat {input.ref_insert_intr} {input.ref_insert_depl} {input.ref_intron_intr} {input.ref_intron_depl} {input.new_ins_intr} {input.new_ins_depl} > output/tmp/all.fa
        awk '/^>/ {{if (seqlen) {{if (seqlen < 1000) print seq; seq=\"\"; seqlen=0}} print}} !/^>/ {{ seq = seq $0; seqlen += length($0) }} END {{if (seqlen < 1000) print seq}}' output/tmp/all.fa > {output.fa}
        touch {output.clusters}
      #  ext/cd-hit-est -i {output.fa} -o {output.clusters} -mask NX -c 0.8 -sc 1 -d 1000 -M 0 -T {threads} -s 0.7
        """
rule filter_fasta:
    input:
        fa=rules.cluster_sequences.output.fa,
    output:
        fa_filtered="output/anc-clust/filtered.fa"
    shell:
        """
        awk 'BEGIN {{ FS="\\n"; RS=">"; ORS=""; count=0; }} 
        NF>1 {{
            header=$1; seq=toupper($2);
            gsub(/\\n/, "", seq);  # Remove newlines within sequences
            if (length(seq) > 0 && length(seq) < 1000) {{
                if (!seqs[header] && !dup_seqs[header,seq]) {{
                    seqs[header] = seq;  # Store the first occurrence of the sequence with this header
                    print ">"header"\\n"seq"\\n";
                }} else if (seqs[header] && seqs[header] != seq) {{
                    # If header is found with a different sequence
                    count = ++dup_counts[header];
                    new_header = header "_alt_haplotype_" count;
                    if (!dup_seqs[header,seq]) {{
                        print ">"new_header"\\n"seq"\\n";  # Print with modified header
                        dup_seqs[header,seq] = 1;  # Mark this sequence as printed
                    }}
                }}
            }}
        }}' {input.fa} > {output.fa_filtered}
        """

rule ancestral_clust:
    threads: 10
    input:
        fa=rules.filter_fasta.output.fa_filtered
    output:
        clust="output/anc-clust/out.clusters"
    shell:
        """
        ext/AncestralClust/ancestralclust -i {input.fa} -o {output.clust} -q output/anc-clust/root.seqs -b 10 -r 100 -p 10 -c 80
        """
rule haplotype_caller:
    input:
        ref_fasta=lambda wildcards: deplete if wildcards.reference == "deplete" else intronerized,
        new_ins = rules.local_denovo.output.vcf,
    output:
        # this vcf has the gatk variants called for each strain, with overlapping
        # records overwritten by denovo assembled insertions (new_ins)
        vcf="output/4-strain-genotyping/{strain}.{reference}.vcf.gz"
    threads: 100
    shell:
        """
        # if [ ! -f {input.ref_fasta}.dict ]; then
        #     gatk CreateSequenceDictionary -R {input.ref_fasta}
        # fi

        samtools view -@ {threads} -h data/strains-bams/{wildcards.strain}.{wildcards.reference}.bam |
            awk 'BEGIN{{OFS="\\t"}} {{if($$0 ~ /^@/) {{print $$0}} else {{n=split($$0,a,"\\t"); s=""; for(i=1;i<=n;i++) {{if(a[i] !~ /^TP:/) {{s=s a[i] OFS}}}} print substr(s, 1, length(s)-1)}}}}' |
            samtools view -@ {threads} -b -o data/strains-bams/{wildcards.strain}.{wildcards.reference}.noTP.bam -
        samtools sort -@ {threads} -o data/strains-bams/{wildcards.strain}.{wildcards.reference}.noTP.sorted.bam data/strains-bams/{wildcards.strain}.{wildcards.reference}.noTP.bam
        samtools index -@ {threads} data/strains-bams/{wildcards.strain}.{wildcards.reference}.noTP.sorted.bam

        picard AddOrReplaceReadGroups \
            I=data/strains-bams/{wildcards.strain}.{wildcards.reference}.noTP.sorted.bam \
            O=data/strains-bams/{wildcards.strain}.{wildcards.reference}.rg.bam \
            RGID=foo \
            RGLB=bar \
            RGPL=illumina \
            RGPU=unit1 \
            RGSM={wildcards.strain} \
            SORT_ORDER=coordinate \
            CREATE_INDEX=True

        ext/gatk-4.5.0.0/gatk MarkDuplicates \
            -I data/strains-bams/{wildcards.strain}.{wildcards.reference}.rg.bam \
            -O data/strains-bams/{wildcards.strain}.{wildcards.reference}.dedup_reads.bam \
            -M data/strains-bams/{wildcards.strain}.{wildcards.reference}.dedup_metrics.txt

        samtools sort -@ {threads} -o data/strains-bams/{wildcards.strain}.{wildcards.reference}.dedup_reads.sorted.bam data/strains-bams/{wildcards.strain}.{wildcards.reference}.dedup_reads.bam
        samtools index -@ {threads} /scratch1/alex/mpop/condon-kramer-2024/data/strains-bams/{wildcards.strain}.{wildcards.reference}.dedup_reads.sorted.bam

        ext/gatk-4.5.0.0/gatk HaplotypeCaller \
            -R {input.ref_fasta} \
            -I data/strains-bams/{wildcards.strain}.{wildcards.reference}.dedup_reads.sorted.bam \
            -O data/tmp/{wildcards.strain}.{wildcards.reference}.tmp.g.vcf.gz

        
        mkdir -p data/tmp/gatk/{wildcards.strain}.{wildcards.reference}
        bcftools index -f data/tmp/{wildcards.strain}.{wildcards.reference}.tmp.g.vcf.gz
        bcftools index -f {input.new_ins}
        bcftools isec -p data/tmp/gatk/{wildcards.strain}.{wildcards.reference} -w 1,2 -O z data/tmp/{wildcards.strain}.{wildcards.reference}.tmp.g.vcf.gz {input.new_ins}
        bcftools concat -a -O z -o data/tmp/{wildcards.strain}.{wildcards.reference}.combined.g.vcf.gz data/tmp/gatk/{wildcards.strain}.{wildcards.reference}/0000.vcf.gz data/tmp/gatk/{wildcards.strain}.{wildcards.reference}/0001.vcf.gz
        bcftools sort -O z -o {output.vcf} data/tmp/{wildcards.strain}.{wildcards.reference}.combined.g.vcf.gz
        bcftools index {output.vcf}
        """

rule filter_variants:
    input:
        vcf_gatk=rules.haplotype_caller.output.vcf,
        ref_fasta=lambda wildcards: deplete if wildcards.reference == "deplete" else intronerized
    output:
        vcf_filtered="output/4-strain-genotyping/{strain}.{reference}.filtered.vcf.gz"
    shell:
        """
        ext/gatk-4.5.0.0/gatk IndexFeatureFile -I {input.vcf_gatk}
        ext/gatk-4.5.0.0/gatk --java-options '-Xmx400g' VariantFiltration -R {input.ref_fasta} -V {input.vcf_gatk} -O {output.vcf_filtered} \
            --filter-expression "QUAL < 30.0" --filter-name "QUAL30"
        """

rule apply_variants:
    input:
        ref=lambda wildcards: deplete if wildcards.reference == "deplete" else intronerized,
        vcf=rules.filter_variants.output.vcf_filtered
    output:
        fasta="output/6-strain-haplotypes/{strain}.{reference}.fasta"
    shell:
        """
        ext/gatk-4.5.0.0/gatk FastaAlternateReferenceMaker -R {input.ref} -V {input.vcf} -O {output.fasta}
        """

rule collect_fastas:
    input:
        fasta=expand("output/6-strain-haplotypes/{strain}.{reference}.fasta", strain=names, reference=["deplete", "intronerized"])
    output:
        ids="output/6-strain-haplotypes/ids.txt"
    shell:
        """
        rm -f {output.ids}
        echo "CCMP1545\t{intronerized}" >> {output.ids}
        echo "RCC1749\t{deplete}" >> {output.ids}
        for f in {input.fasta}; do
            strain=$(basename $f | cut -d. -f1)
            ref=$(basename $f | cut -d. -f2)
            id="{$strain}_v_{$ref}"
            echo "$id\t$f" >> {output.ids}
        done
        """

rule blast_all_vs_all:
    input:
        ref_introns = lambda wildcards: known_introns_deplete if wildcards.reference == "deplete" else known_introns_intronerized
    output: "output/blast.knownintrons.{reference}.tsv"
    shell:
        """
        makeblastdb -in {input.ref_introns} -dbtype nucl
        blastn -num_threads 100 -query {input.ref_introns} -db {input.ref_introns} -outfmt '6 qseqid sseqid pident length qlen slen' -out {output} -max_target_seqs 20
        """

rule no_overlap:
    input: lambda wildcards: known_introns_deplete if wildcards.reference == "deplete" else known_introns_intronerized
    output: "output/pre.knownintrons.{reference}.fa"
    run:
        contig_map = defaultdict(list)
        for record in SeqIO.parse(input[0], 'fasta'):
            contig = record.id.split('.')[3].split(':')[0]
            start = int(record.id.split(':')[1].split('-')[0])
            stop = int(record.id.split(':')[1].split('-')[1])
            contig_map[contig].append((start, stop, record))
        
        longest_records = []

        # For each contig, find non-overlapping sequences
        for contig, intervals in contig_map.items():
            intervals.sort()  # Sort intervals by start position
            non_overlapping = []

            current_interval = intervals[0]
            for interval in intervals[1:]:
                if interval[0] > current_interval[1]:  # No overlap
                    non_overlapping.append(current_interval[2])
                    current_interval = interval
                else:  # Overlap, choose the longer sequence
                    if len(interval[2]) > len(current_interval[2]):
                        current_interval = interval
            
            non_overlapping.append(current_interval[2])
            longest_records.extend(non_overlapping)
        with open(output[0], 'w') as out_f:
            for record in longest_records:
                SeqIO.write(record, out_f, 'fasta')

rule filter_hits:
    input: rules.blast_all_vs_all.output, rules.no_overlap.output
    output: "output/filtered.introners.from_intron.{reference}.fa"
    run:
        import pandas as pd
        from Bio import SeqIO

        # Load BLAST results
        df = pd.read_csv(input[0], sep='\t', header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen'])
        
        # Filter for hits that meet criteria
        df = df[(df['qseqid'] != df['sseqid']) & (df['pident'] > 90) & (df['length'] / df['qlen'] > 0.9) & (df['length'] / df['slen'] > 0.9) & (df['length'] < 10000)]

        # Get unique sequences that have valid hits
        valid_sequences = set(df['qseqid'].unique())


        # Read input sequences and write only the valid ones to the output
        with open(output[0], 'w') as out_f:
            for record in SeqIO.parse(input[1], 'fasta'):
                if record.id in valid_sequences:
                    SeqIO.write(record, out_f, 'fasta')


rule indels_vs_introners:
    input: q = lambda wildcards: rules.prepare_indels.output.depl if wildcards.reference == "deplete" else rules.prepare_indels.output.intr,
              s = rules.filter_hits.output
    output: "output/blast.indels.{reference}.tsv"
    shell:
        """
        makeblastdb -in {input.s} -dbtype nucl
        blastn -num_threads 100 -query {input.q} -db {input.s} -outfmt '6 qseqid sseqid pident length qlen slen' -out {output} -max_target_seqs 20
        """

rule filter_hits_indels_vs_introners:
    input: lambda wildcards: rules.prepare_indels.output.depl if wildcards.reference == "deplete" else rules.prepare_indels.output.intr, rules.indels_vs_introners.output
    output: "output/filtered.introners.from_indels.{reference}.fa"
    run:
        import pandas as pd
        from Bio import SeqIO

        # Load BLAST results
        df = pd.read_csv(input[1], sep='\t', header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen'])
        
        # Filter for hits that meet criteria
        df = df[(df['pident'] > 90) & (df['length'] / df['qlen'] > 0.9) & (df['length'] / df['slen'] > 0.9) & (df['length'] < 10000)]

        # Get unique sequences that have valid hits
        valid_sequences = set(df['qseqid'].unique())

        # Read input sequences and write only the valid ones to the output
        with open(output[0], 'w') as out_f:
            for record in SeqIO.parse(input[0], 'fasta'):
                if record.id in valid_sequences:
                    SeqIO.write(record, out_f, 'fasta')

rule cat_indel_introners_and_intron_introners:
    input:
        indel = rules.filter_hits_indels_vs_introners.output,
        intron = rules.filter_hits.output
    output: "output/filtered.introners.from_intron_and_indel.{reference}.fa"
    shell:
        """
        cat {input.indel} {input.intron} > {output}
        """

rule localdenovo_vs_introners:
    input: q = rules.cat_fa.output.all_fa, s = rules.cat_indel_introners_and_intron_introners.output
    output: "output/blast.denovo_vs_introners.{reference}.tsv"
    shell:
        """
        makeblastdb -in {input.s} -dbtype nucl
        blastn -num_threads 100 -query {input.q} -db {input.s} -outfmt '6 qseqid sseqid pident length qlen slen' -out {output} -max_target_seqs 20
        """


rule filter_hits_denovo_vs_introners:
    input: rules.localdenovo_vs_introners.output, rules.cat_fa.output.all_fa
    output: "output/filtered.introners.from_denovo.{reference}.fa"
    run:
        import pandas as pd
        from Bio import SeqIO

        # Load BLAST results
        df = pd.read_csv(input[0], sep='\t', header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen'])
        
        # Filter for hits that meet criteria
        df = df[(df['pident'] > 90) & (df['length'] / df['qlen'] > 0.9) & (df['length'] / df['slen'] > 0.9) & (df['length'] < 10000)]

        # Get unique sequences that have valid hits
        valid_sequences = set(df['qseqid'].unique())

        # Read input sequences and write only the valid ones to the output
        with open(output[0], 'w') as out_f:
            for record in SeqIO.parse(input[1], 'fasta'):
                if record.id in valid_sequences:
                    SeqIO.write(record, out_f, 'fasta')

rule cat_indel_introners_and_intron_introners_and_denovo_introners:
    input:
        indels_and_intron=rules.cat_indel_introners_and_intron_introners.output,
        denovo=rules.filter_hits_denovo_vs_introners.output
    output: "output/filtered.introners.from_intron_and_indel_and_denovo.{reference}.fa"
    shell:
        """
        cat {input.indels_and_intron} {input.denovo} > {output}
        """

rule prepare_graph_indels:
    input:
        indels = "mega/mpusilla_mega_graph.vcf"
    output:
        depl="mega/insertions-in-deplete.fa",
        intr="mega/insertions-in-intronerized.fa"
    run:
        with open(input.indels, 'r') as f, \
            open(output.depl, 'w') as depl, open(output.intr, 'w') as intr:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                scaf, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
                if len(ref) > len(alt) and len(ref) > 40:
                    intr.write(f">MegaGraph.Insertion.Intronerized.{scaf}:{pos}\n{ref}\n")
                elif len(alt) > len(ref) and len(alt) > 40:
                    depl.write(f">MegaGraph.Insertion.Deplete.{scaf}:{pos}\n{alt}\n")



rule graph_indels_vs_introners:
    input: q = lambda wildcards: rules.prepare_graph_indels.output.depl if wildcards.reference == "deplete" else rules.prepare_graph_indels.output.intr,
              s = rules.cat_indel_introners_and_intron_introners_and_denovo_introners.output
    output: "output/blast.graph_indels.{reference}.tsv"
    shell:
        """
        makeblastdb -in {input.s} -dbtype nucl
        blastn -num_threads 100 -query {input.q} -db {input.s} -outfmt '6 qseqid sseqid pident length qlen slen' -out {output} -max_target_seqs 20
        """

rule filter_hits_graph_indels_vs_introners:
    input: lambda wildcards: rules.prepare_graph_indels.output.depl if wildcards.reference == "deplete" else rules.prepare_graph_indels.output.intr, rules.graph_indels_vs_introners.output
    output: "output/filtered.introners.from_megagraph_indels.{reference}.fa"
    run:
        import pandas as pd
        from Bio import SeqIO

        # Load BLAST results
        df = pd.read_csv(input[1], sep='\t', header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen'])
        
        # Filter for hits that meet criteria
        df = df[(df['pident'] > 90) & (df['length'] / df['qlen'] > 0.9) & (df['length'] / df['slen'] > 0.9) & (df['length'] < 10000)]

        # Get unique sequences that have valid hits
        valid_sequences = set(df['qseqid'].unique())

        # Read input sequences and write only the valid ones to the output
        with open(output[0], 'w') as out_f:
            for record in SeqIO.parse(input[0], 'fasta'):
                if record.id in valid_sequences:
                    SeqIO.write(record, out_f, 'fasta') 
rule cleanup:
    input: expand('output/filtered.introners.from_{source}.{reference}.fa', source=['intron', 'indels', 'denovo', 'megagraph_indels'], reference=['deplete', 'intronerized'])
    output: 'output/calls/filtered.introners.from_{source}.{reference}.fa'
    shell:
        """
        cat "output/filtered.introners.from_{wildcards.source}.{wildcards.reference}.fa" > {output}
        """