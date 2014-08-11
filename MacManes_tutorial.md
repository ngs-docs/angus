Assembly with SOAPdenovo-Trans
--

Startup AMI

	sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
        default-jre pkg-config libncurses5-dev r-base-core \
        r-cran-gplots python-matplotlib sysstat samtools python-biopython

Installation

	wget http://downloads.sourceforge.net/project/soapdenovotrans/SOAPdenovo-Trans/bin/v1.03/SOAPdenovo-Trans-bin-v1.03.tar.gz

	tar -zxf SOAPdenovo-Trans-bin-v1.03.tar.gz
		

make SOAP folder and get into it
	
	mkdir SOAP
	cd SOAP
	
Make config file
	
	nano config.txt

    #maximal read length
    max_rd_len=100
    [LIB]
    #maximal read length in this lib
    rd_len_cutof=45
    #average insert size
    avg_ins=200
    #if sequence needs to be reversed
    reverse_seq=0
    #in which part(s) the reads are used
    asm_flags=3
   	#minimum aligned length to contigs for a reliable read location 
    map_len=32
    #fastq file for read 1
    q1=/path/**LIBNAMEA**/fastq_read_1.fq
    #fastq file for read 2 always follows fastq file for read 1
    q2=/path/**LIBNAMEA**/fastq_read_2.fq
    #fasta file for read 1
    q=/path/**LIBNAMEA**/fastq_read_single.fq
 


Assembly optimization
	
	mkdir SOAP
	nano config.txt
	for k in 31 41 51 61 71 91;
		do SOAPdenovo-Trans-127mer all -L 300 -p 4 -K $k -s config.txt -o assembly$k; done
		
	
Pick the best assembly

	Transrate -> http://hibberdlab.com/transrate/


Install and run cd-hit est
	
	wget https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
	tar -zxf cd-hit-v4.6.1-2012-08-27.tgz
	cd cd-hit-v4.6.1-2012-08-27
	make
	PATH=$PATH:home/ubuntu/cd-hit-v4.6.1-2012-08-27
	cd-hit-est -i Trinity_all_X.fasta -o trin.fasta

bwa -> 

	git clone https://github.com/lh3/bwa.git
	cd bwa
	make
	PATH=$PATH:/home/ubuntu/bwa


eXpress -> http://bio.math.berkeley.edu/eXpress/

	wget http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz

	tar -xzf express-1.5.1-linux_x86_64.tgz
	
	cd express-1.5.1-linux_x86_64/
	
	PATH=$PATH:/home/ubuntu/express-1.5.1-linux_x86_64
	

Do Mapping


	mkdir soap_index && cd soap_index
	
	bwa index -p all ../../soap_assemblies/assembly51.fasta
	
    bwa mem -t 4 soap_all_index/all \
    ../trimmed_x/ORE_sdE3_rep1_1_pe \
    ../trimmed_x/ORE_sdE3_rep1_2_pe  | \
    express -o express_soap/ORE_sdE3_rep1.xprs \
    -p4 soap_all_index/soap.fasta

