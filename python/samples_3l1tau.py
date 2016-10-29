from collections import OrderedDict as OD

samples = OD()

samples["ttHJetToNonbb_M125"] = OD([
  ("type",                  "mc"),
  ("sample_category",       "signal"),
  ("process_name",          "ttHToNonbb_M125"),
  ("process_name_specific", OD([
    ("2015", ["ttHToNonbb_M125", ]),
    ("2016", ["ttHToNonbb_M125", ]),
  ])),
  ("use_it",                True),
  ("xsection",              0.2151),
])

samples["TTZToLLNuNu"] = OD([
  ("type",                  "mc"),
  ("sample_category",       "TTZ"),
  ("process_name",          "TTZToLLNuNu"),
  ("process_name_specific", OD([
    ("2015", ["TTZToLLNuNu", ]),
    ("2016", ["TTZToLLNuNu_M-10",
              "TTZToLLNuNu_M-10_ext1", ]),
  ])),
  ("use_it",                True),
  ("xsection",              0.2529),
])
