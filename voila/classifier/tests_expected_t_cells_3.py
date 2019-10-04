sg_file = '/home/paul/PycharmProjects/majiq/test_cases/classifier/t_cells_2/splicegraph.sql'
psi_file = '/home/paul/PycharmProjects/majiq/test_cases/classifier/t_cells_2/ran_treg.psi.voila'
sg_file = '/Users/calebradens/Documents/majiq_dev/classifier_dev/splicegraph.sql'
psi_file = '/Users/calebradens/Documents/majiq_dev/classifier_dev/ran_treg.psi.voila'
# source ~/bin/majiq2/voila_viewer/voila_view_env/bin/activate
# voila view /Users/calebradens/Documents/majiq_dev/classifier_dev/splicegraph.sql /Users/calebradens/Documents/majiq_dev/classifier_dev/ran_treg.psi.voila

expected_modules = {
#"gene:ENSG00000083168": [
    #     {
    #         "cassette_exon": "1",
    #         "alt3ss": "",
    #         "alt5ss": "1",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "",
    #         "complex": "True"
    #     },
    #     {
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         # "p_alt5ss": "1",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "1",
    #         "alternative_intron": "2",
    #         "complex": "True"},
    #     {
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "1",
    #         "p_alt3ss": "",
    #         # "p_alt5ss": "1",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "",
    #         "complex": "False"
    #     },
    #     {},
    #     {
    #         "cassette_exon": "1",
    #         "alt3ss": "1",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "",
    #         "complex": "True"
    #     }
    # ],
# "gene:ENSG00000081307": [
    #     {
    #         "module_id": "gene:ENSG00000081307_1",
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "1",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "1",
    #         "alternative_intron": "2",
    #         "complex": "True"
    #     },
    #     {
    #         "module_id": "gene:ENSG00000081307_2",
    #         "cassette_exon": "1",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "1",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "1",
    #         "complex": "True"
    #     }
    # ],
# "gene:ENSG00000111653": [
#         {
#             "module_id": "gene:ENSG00000111653_1",
#             "cassette_exon": "2",
#             "alt3ss": "",
#             "alt5ss": "",
#             "p_alt3ss": "",
#             "p_alt5ss": "",
#
#             "mutually_exclusive": "",
#             "ale": "",
#             "afe": "",
#             "p_ale": "",
#             "p_afe": "",
#             "multi_exon_spanning": "2",
#             "alternative_intron": "",
#             "complex": "True"
#         },
#         {
#             "module_id": "gene:ENSG00000111653_2",
#             "cassette_exon": "",
#             "alt3ss": "",
#             "alt5ss": "1",
#             "p_alt3ss": "",
#             "p_alt5ss": "",
#             "mutually_exclusive": "",
#             "ale": "",
#             "afe": "",
#             "p_ale": "",
#             "p_afe": "",
#             "multi_exon_spanning": "",
#             "alternative_intron": "",
#             "complex": "False"
#         },
#         {
#             "module_id": "gene:ENSG00000111653_3",
#             "cassette_exon": "",
#             "alt3ss": "",
#             "alt5ss": "",
#             "p_alt3ss": "",
#             "p_alt5ss": "",
#
#             "mutually_exclusive": "",
#             "ale": "1",
#             "afe": "",
#             "p_ale": "",
#             "p_afe": "",
#             "multi_exon_spanning": "",
#             "alternative_intron": "",
#             "complex": "False"
#         }
#     ],
"gene:ENSG00000082074": [

        {
            "module_id": "gene:ENSG00000082074_1",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000082074_2",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "1",
            "p_alt5ss": "1",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "2",
            "complex": "True"
        },
        {
            "module_id": "gene:ENSG00000082074_3",
            "alternative_intron": "1",
        },
        {
            "module_id": "gene:ENSG00000082074_4",
            "cassette_exon": "1",
            "alternative_intron": "2",
            "alt5ss": "1",
            "p_alt3ss": "1",
            "p_alt5ss": "2",
            "p_ale": "1",
        },
        {
            "module_id": "gene:ENSG00000082074_5",
            "alternative_intron": "1",
        },
        {
            "module_id": "gene:ENSG00000082074_6",
            "alternative_intron": "1",
        },
        {
            "module_id": "gene:ENSG00000082074_7",
            "alternative_intron": "1",
        },
        {
            "module_id": "gene:ENSG00000082074_8",
            "alternative_intron": "1",
        },
    ],
    "gene:ENSG00000081721": [
        {
            "module_id": "gene:ENSG00000081721_1",
            "lsv_ids": "gene:ENSG00000081721:s:161749758-161750145;gene:ENSG00000081721:t:161751655-161751781",
            "cassette_exon": "",
            "alt3ss": "1",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        }
    ],

    "gene:ENSG00000177239": [
        {"cassette_exon": "1",
         "alt3ss": "",
         "alt5ss": "",
         "p_alt3ss": "",
         "p_alt5ss": "",
         "mutually_exclusive": "",
         "ale": "",
         "afe": "",
         "p_ale": "",
         "p_afe": "",
         "multi_exon_spanning": "",
         "alternative_intron": "",
         },
        {"cassette_exon": "",
         "alt3ss": "",
         "alt5ss": "",
         "p_alt3ss": "",
         "p_alt5ss": "",
         "mutually_exclusive": "",
         "ale": "",
         "afe": "",
         "p_ale": "",
         "p_afe": "",
         "multi_exon_spanning": "",
         "alternative_intron": "1",
         },
        {"cassette_exon": "1",
         "alt3ss": "1",
         "alt5ss": "",
         "alt3and5ss": "2",
         "p_alt3ss": "",
         "p_alt5ss": "",
         "mutually_exclusive": "",
         "ale": "1",
         "afe": "",
         "p_ale": "",
         "p_afe": "",
         "multi_exon_spanning": "1",
         "alternative_intron": "1",
         "exitron": ""  #hmmph
         },
        {"cassette_exon": "",
         "alt3ss": "1",
         "alt5ss": "",
         "p_alt3ss": "",
         "p_alt5ss": "",
         "mutually_exclusive": "",
         "ale": "",
         "afe": "",
         "p_ale": "",
         "p_afe": "",
         "multi_exon_spanning": "",
         "alternative_intron": "",
         },

    ],
    "gene:ENSG00000205744": [
        {}, {}, {},
        {
            "mutually_exclusive": "1",
            "alternative_intron": ""
        }
    ],
    "gene:ENSG00000003756": [
        {}, {}, {}, {}, {}, {}, {}, {}, {},
        {
            'tandem_cassette': '1'
        },
        {}, {}
    ],
    "gene:ENSG00000164674":[
        {},
        {
            "module_id": "gene:ENSG00000164674_2",
            "mutually_exclusive": "",
            # "exitron": "1",
            "p_alt3ss": "",
            "p_alt5ss": "",
            "alternative_intron": "1"

        }, {}, {}, {}, {}
    ],
    "gene:ENSG00000111540": [
        {
            "module_id": "gene:ENSG00000111540_1",
            "lsv_ids": "gene:ENSG00000111540:s:55974011-55974139",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "2",
            "afe": "1",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "1",
            "alternative_intron": "",
            "complex": "True"
        }, {
            "ale": "1",
            "afe": "1"
        }
    ],
    "gene:ENSG00000266967": [
        {
            'afe': '1'
        }, {}
    ],
    "gene:ENSG00000134352": [
        {
            'tandem_cassette': '1'
        }, {}, {}, {}, {}, {}, {}
    ],
    "gene:ENSG00000170581": [
        {}, {}, {}, {}, {
            'tandem_cassette': ''
        }, {}
    ],
    "gene:ENSG00000082213": [
        {
            "module_id": "gene:ENSG00000082213_1",
            "lsv_ids": "gene:ENSG00000082213:s:31535744-31535893",
            "cassette_exon": "",
            "alt3ss": "1",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000082213_2",
            "lsv_ids": "gene:ENSG00000082213:s:31545646-31545714;gene:ENSG00000082213:t:31551293-31551432",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        }
    ],
    "gene:ENSG00000047249": [
        {
            "module_id": "gene:ENSG00000047249_1",
            "lsv_ids": "gene:ENSG00000047249:t:53801799-53801896;gene:ENSG00000047249:s:53811164-53811217",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000047249_2",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            # "ale": "",
            # "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "1",
            "alternative_intron": "",
            # "complex": "False"
        }
    ],
    "gene:ENSG00000082898": [
        {
            "module_id": "gene:ENSG00000082898_1",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            # "afe": "1",
            "p_ale": "",
            "p_afe": "",
            # "multi_exon_spanning": "1",
            # "alternative_intron": "1",
            # "complex": "True"
        }
    ],

    "gene:ENSG00000175582": [
        {
            "module_id": "gene:ENSG00000175582_1",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "1",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        }
    ],
    "gene:ENSG00000144535": [
        {
            "module_id": "gene:ENSG00000144535_1",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "1",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000144535_2",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000144535_3",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        }
    ],

    # "gene:ENSG00000285708": [
    #     {
    #         "module_id": "gene:ENSG00000285708_1",
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "1",
    #         "complex": "False"
    #     },
    #     {
    #         "module_id": "gene:ENSG00000285708_2",
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "1",
    #         "complex": "False"
    #     },
    #     {
    #         "module_id": "gene:ENSG00000285708_3",
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "1",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "",
    #         "complex": "False"
    #     },
    #     {
    #         "module_id": "gene:ENSG00000285708_4",
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "1",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "",
    #         "complex": "False"
    #     },
    #     {
    #         "module_id": "gene:ENSG00000285708_5",
    #         "cassette_exon": "",
    #         "alt3ss": "1",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "",
    #         "complex": "False"
    #     },
    #     {
    #         "module_id": "gene:ENSG00000285708_6",
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "1",
    #         "complex": "False"
    #     },
    #     {
    #         "module_id": "gene:ENSG00000285708_7",
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "1",
    #         "complex": "False"
    #     },
    #     {
    #         "module_id": "gene:ENSG00000285708_8",
    #         "cassette_exon": "",
    #         "alt3ss": "",
    #         "alt5ss": "",
    #         "p_alt3ss": "",
    #         "p_alt5ss": "",
    #
    #         "mutually_exclusive": "",
    #         "ale": "",
    #         "afe": "",
    #         "p_ale": "",
    #         "p_afe": "",
    #         "multi_exon_spanning": "",
    #         "alternative_intron": "1",
    #         "complex": "False"
    #     }
    # ],
    "gene:ENSG00000151327": [
        {
            "module_id": "gene:ENSG00000151327_1",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "1",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "afe": "1",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "True"
        },
        {
            "module_id": "gene:ENSG00000151327_2",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "1",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "afe": "1",
            "ale": "1",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "1",
            "alternative_intron": "1",
            "complex": "True"
        },
        {
            "module_id": "gene:ENSG00000151327_3",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "1",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        }
    ],
    "gene:ENSG00000144381": [
        {
            "module_id": "gene:ENSG00000144381_1",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "2",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "1",
            "alternative_intron": "",
            "complex": "True"
        },
        {
            "module_id": "gene:ENSG00000144381_2",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "1",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000144381_3",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "1",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000144381_4",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "1",
            "p_afe": "1",
            "multi_exon_spanning": "",
            "alternative_intron": "1",
            "complex": "True"
        },
        {
            "module_id": "gene:ENSG00000144381_5",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "1",
            "complex": "False"
        }
    ],
    "gene:ENSG00000090863": [
        {
            "module_id": "gene:ENSG00000090863_1",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000090863_2",
            "cassette_exon": "2",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "1",
            "p_alt5ss": "1",

            "mutually_exclusive": "1",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "1",
            "alternative_intron": "1",
            "complex": "True"
        },
        {
            "module_id": "gene:ENSG00000090863_3",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "1",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000090863_4",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "1",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000090863_5",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "1",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "1",
            "alternative_intron": "2",
            "complex": "True"
        },
        {
            "module_id": "gene:ENSG00000090863_6",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "1",
            "complex": "False"
        }
    ],
    "gene:ENSG00000160201": [
        {
            "module_id": "gene:ENSG00000160201_1",
            "cassette_exon": "4",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "1",
            "alternative_intron": "",
            "complex": "True"
        },
        {
            "module_id": "gene:ENSG00000160201_2",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000160201_3",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "1",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000160201_4",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "1",
            "p_alt5ss": "1",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "2",
            "complex": "True"
        }
    ],

    "gene:ENSG00000081665": [
        {
            "module_id": "gene:ENSG00000081665_1",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "1",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000081665_2",
            "cassette_exon": "",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "1",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "1",
            "complex": "True"
        }
    ],
    "gene:ENSG00000111412": [
        {
            "module_id": "gene:ENSG00000111412_1",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "1",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "1",
            "alternative_intron": "",
            "complex": "True"
        }
    ],
    "gene:ENSG00000162735": [
        {
            "module_id": "gene:ENSG00000162735_1",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        },
        {
            "module_id": "gene:ENSG00000162735_2",
            "cassette_exon": "1",
            "alt3ss": "",
            "alt5ss": "",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "False"
        }
    ],
    "gene:ENSG00000119397": [
        {}, {}, {}, {}, {},
        {
            "module_id": "gene:ENSG00000119397_6",
            "tandem_cassette": "1",
            "multi_exon_spanning": "1",
            "complex": "False"
        }, {}, {}, {}
    ],
    'gene:ENSG00000114388': [
        {
            "module_id": "gene:ENSG00000114388_1",
            "cassette_exon": "",
            "alt3ss": "1",
            "alt5ss": "1",
            "alt3and5ss": "1",
            "p_alt3ss": "",
            "p_alt5ss": "",

            "mutually_exclusive": "",
            "ale": "",
            "afe": "",
            "p_ale": "",
            "p_afe": "",
            "multi_exon_spanning": "",
            "alternative_intron": "",
            "complex": "True"
        },
        {}
    ],
    'gene:ENSG00000236449': [
        {
            "module_id": "gene:ENSG00000236449_1",
            "orphan_junction": "1",
            "p_ale": "2",
            "p_afe": "",
        }
    ],
    'gene:ENSG00000233175': [
        {
            "module_id": "gene:ENSG00000233175_1",
            # "exitron": "2"
        },

    ],
    "gene:ENSG00000122390":[
        {}, {}, {}, {}, {}
    ],
    "gene:ENSG00000087206":[
        {},{},{}
    ],
    "gene:ENSG00000100796":[

    ],
    "gene:ENSG00000124151":[

    ]

}

expected_modules_constitutive = {
    'gene:ENSG00000218739': [
        {}, {}, {}, {
            'constitutive_intron': '1'
        }

    ],
    'gene:ENSG00000273189': [
        {
            'constitutive_junction': '1'
        }, {
            'constitutive_junction': '1'
        }, {
            'constitutive_junction': '1'
        }, {
            'constitutive_junction': '1'
        }, {
            'constitutive_junction': '1'
        }, {
            'constitutive_junction': '1'
        }, {
            'constitutive_intron': '1'
        }
    ],
    "gene:ENSG00000223509":[
        {},
        {'constitutive_junction': '1'},
        {'constitutive_junction': '1'},
        {'constitutive_junction': '1'},
        {},
        {}
    ]
}

expected_cassette_exons = {
    'gene:ENSG00000138073': [
        {
            'ran_treg_E(PSI)': '0.125',
        },{
            'ran_treg_E(PSI)': '0.875',
        },{
            'ran_treg_E(PSI)': '0.07913184',
        },{
            'ran_treg_E(PSI)': '0.90108526',
        }
    ]
}

expected_alternative_intron = {
    'gene:ENSG00000123146': [
        {'Reference Exon Coordinate': '14381156-14381545'},
        {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
    ]
}

expected_alt3ss = {
    "gene:ENSG00000100796":[
        {'ran_treg_E(PSI)': '0.992'},
        {'ran_treg_E(PSI)': '0.008'}
    ],
    "gene:ENSG00000124151":[
        {'ran_treg_E(PSI)': '0.86'},
        {'ran_treg_E(PSI)': '0.14'},
        {'ran_treg_E(PSI)': '0.665'},
        {'ran_treg_E(PSI)': '0.335'}
    ]
}

expected_alt5ss = {
    "gene:ENSG00000122390":[
        {'ran_treg_E(PSI)': '0.391'},
        {'ran_treg_E(PSI)': '0.609'}
    ],
    "gene:ENSG00000087206":[
        {'ran_treg_E(PSI)': '0.039'},
        {'ran_treg_E(PSI)': '0.961'}
    ]
}

# expected_alt5ss = {
#     'gene:ENSG00000138085'
# }

expected_mpes = {
    "gene:ENSG00000218739":[
        {"LSV ID(s)":"",
         "Reference Exon Constant Region":"37196488-37196520"},
        {"LSV ID(s)":"",
         "Reference Exon Constant Region":"37199704-37199819"},
        {"LSV ID(s)":"gene:ENSG00000218739:t:37201048-37201092",
         "Reference Exon Constant Region":"37201048-37201092"},
        {"LSV ID(s)":"gene:ENSG00000218739:s:37201642-37201726",
         "Reference Exon Constant Region":"37201642-37201726"},
        {"LSV ID(s)":"gene:ENSG00000218739:t:37204304-37204741",
         "Reference Exon Constant Region":"37204304-37204741"}
    ],
    "gene:ENSG00000223509":[
        {"Constitutive Exon or Intron":"e;e;e"},
        {"Constitutive Exon or Intron":"h"},#half exon
        {"Constitutive Exon or Intron":"e;e;e",
         "Constitutive De Novo":"False;False;False"}
    ]
}