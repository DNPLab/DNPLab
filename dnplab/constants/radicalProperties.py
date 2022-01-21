#######################################
# EPR Properties of Selected Radicals #
#######################################
#
radicalProperties = {}

# When adding new radicals, make sure to add a literature reference for each entry. Keys in dictonary should all be lower case.

# Free electron
# http://physics.nist.gov/constants
radicalProperties["gfree"] = [2.00231930436153, None, [0]]


# TEMPO in Toluene
# Table 1, page 243 in Lebedev et al., Bioactive Spinlabels, Springer Verlag, 1992
radicalProperties["tempo1"] = [[2.00980, 2.00622, 2.00220], "14N", [16.8, 20.5, 95.9]]

# TEMPO in Methanol
# Table 1, page 243 in Lebedev et al., Bioactive Spinlabels, Springer Verlag, 1992
radicalProperties["tempo2"] = [[2.00909, 2.00621, 2.00222], "14N", [20.2, 20.2, 102.1]]

# BDPA in Polystyrene
# Bennati et al., JMR, 1999. (2H couplings were scaled to 1H)
radicalProperties["bdpa"] = [[2.00263, 2.00260, 2.00257], "1H", [50.2, 34.5, 13.0]]
