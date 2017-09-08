def isPairedSRA(filename):
    filename = os.path.abspath(filename);
    try:
        contents = sp.check_output(["/home/yunjie/Downloads/Software/sratoolkit.2.6.2-ubuntu64/bin/fastq-dump","-X","1","-Z","--split-spot", filename]);
    except sp.CalledProcessError, e:
        raise Exception("Error running fastq-dump on",filename);

    if(contents.count("\n") == 4):
        return False;
    elif(contents.count("\n") == 8):
        return True;
    else:
        raise Exception("Unexpected output from fast-dump on ", filename);
