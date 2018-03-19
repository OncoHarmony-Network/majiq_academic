void func1 (string bamfile){



    SeqLib::BamReader br;
    br.Open(bamfile);


    SeqLib::GRC g(GZBED, br.Header());

    BOOST_CHECK_EQUAL(g.size(), 3);

    BOOST_CHECK_EQUAL(g[2].chr, 21);

    SeqLib::GRC v(GZVCF, br.Header());
    BOOST_CHECK_EQUAL(v.size(), 57);

    BOOST_CHECK_EQUAL(v[29].chr, 0);

}