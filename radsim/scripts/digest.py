
def main(opts):
    """Does the actual digestions given an option namespace"""

    # gets enzyme class by name
    cutters = opts.enzyme.split(",")
    cutters = [getattr(Restriction, enz) for enz in cutters]
    print("Using the following exzymes:")
    for cutter in cutters:
        print("\t{}".format(cutter))
        print("\t{}".format(cutter.site))

    # Opens file, and creates fasta reader
    seq_file = open(opts.seqfile, "rb")
    seqs = SeqIO.parse(seq_file, "fasta")

    # Digest all sequences in the fasta file
    count = 0
    for record in seqs:
        record_count = 0
        # When we're counting, we only want to show how many cut sites there
        # are. Given that cutting a sequence 0 times gives one fragment, we
        # need to decrement this, or we add one to the true number of sites

        # Do virtual digest
        # TODO: change this to use `cutters` and RestrictionBatch
        fragments = cutter.catalyse(record.seq)

        # Find fragment lenghts
        fragment_lengths = []
        for seq in fragments:
            fragment_lengths.append(len(seq))

        # Count how many fragments meet the specified selection criteria
        for length in fragment_lengths:
            if opts.count:  # Count everything
                record_count += max(0, len(fragment_lengths) - 1)
            elif length > opts.minlen and  length < opts.maxlen:
                record_count += 1

        # Append the counts for this record to the total sum
        count += record_count

        if opts.verbose > 0:
            # Print summary
            print "%s has %i RADseq tags" % (record.id, record_count)

    print "In total, %i RADseq tags were found" % count
