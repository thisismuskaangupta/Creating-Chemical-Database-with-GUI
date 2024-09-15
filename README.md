# Creating-Chemical-Database-with-GUI
this contains the code written for the MSc Bioinformatics course - Using Chemical Structure Databases in Drug Discovery. grade: A3:20 in associated presentation and viva; A4:19 in associated written report.

this script takes an SDF file (containing information on molecules), creates a database of these molecules (named '2873826.db'), and then uses this database to create a GUI to access and search this database.

names are hardcoded as the focus of this assignment was a working GUI.

Usage - python3 script_name.py

Please ensure that you have three scripts – ‘2873826_DatabasePyScript.py’, ‘2873826_GUIPyScript.py’, and ‘2873826_SQL_DDL_Script.sql’ in the current working directory. Please also ensure you have ‘Molecules8.sdf’ in the current working directory. The scripts have sensible error handling to direct the user on how to use them. First, the Database script must be run, and then the GUI script must be run.

Extra Notes - 
1) the following characteristics are inserted in the database, some extracted directly from the sdf, others calculated using external libraries or self-written functions.
  a) QED,
  b) SPS,
  c) Molecular Weight,
  d) Exact Molecular Weight, 
  e) Number of Valence Electrons,
  f) Number of Radical Electrons,
  g) H Bond Donors,
  h) H Bond Acceptors,
  i) LogP,
  j) LogD,
  k) Number of Aromatic Rings,
  l) Molecular Formula,
  m) Name of Compounds,
  n) Ring Count, 
  o) Rotatable Bond Count,
  p) TPSA,
  q) Molar Refractivity,
  r) Number of Atoms,
  s) Formal Charge,
  t) Number of Heavy Atoms,
  u) Molecular Structure,
  v) SMILES code
  w) fused aromatic ring count, and
  x) all filters (Lipinski, lead-likeness, bioavailability, Ghose, Veber, REOS, drug-like QED, rule of 3). (for each filter, for every metric that was fulfilled, a score of 1 was assigned. E.g. Lipinski’s filter       contains 4 metrics, and if a compound satisfies all 4, then it passes the filter. Instead of using Booleans, a scoring system was used to count how many metrics were satisfied for each filter, as this       contains more information than a true/false binary.)
2)  After desired searching and filtering is done in the GUI, the names of the set of resultant compounds can be exported as a CSV file, named ‘results.csv’, using the ‘Export’ button. 
