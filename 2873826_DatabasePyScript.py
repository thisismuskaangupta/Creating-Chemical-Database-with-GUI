#importing logging
try:
    import logging
except ModuleNotFoundError:
    print("Hello, the logging module is not installed. Please install it and run the script again. Raising SystemExit.")
    raise SystemExit(1)

#setting up the logger
#getting the logger
logger = logging.getLogger('logger')
#setting the level of the logger
logger.setLevel(logging.WARNING)
#creating a stream handler.
shandler = logging.StreamHandler()
#formatting the stream handler
shandler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - {%(pathname)s:%(lineno)d} - %(message)s'))
#adding the stream handler to the logger.
logger.addHandler(shandler)

#importing necessary modules
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import Descriptors
    import sqlite3
    from STOUT import translate_forward
    import os
except ModuleNotFoundError as e:
    logger.error(f"Hello, the following library was not found. Please install it and run the script again, Raising SystemExit.\n{e}")
    raise SystemExit(1)

#setting file paths
database_path = '2873826.db'
sdf_filepath = 'Molecules8.sdf'

#making sure the sdf database filepath exists
try:
    assert os.path.isfile(sdf_filepath) == True
except AssertionError:
    logger.error(f"Hello the given sdf file {sdf_filepath} does not exist. Please fix the error or make sure the file exists. Raising SystemExit.")
    raise SystemExit(1)

#creating an rdkit database of the sdf file
db = Chem.SDMolSupplier(sdf_filepath)

#defining functions 

#defining a function to calculate the number of the fused aromatic rings in a compounds
def fused_aro_rings_calculator(list_supplied): #as an argument, it takes a list of lists where each list is a tuple of the atom numbers in an aromatic ring
    dictionary_results = {} #starting with an empty dictionary to store the results in
    fuses = 0 #starting with a zero number of fuses counted, each fuse is defined as a common bonds between any two aromatic rings in the input list supplied.
    for item in list_supplied:
        dictionary_results[str(item)] = 0 #using the ring as a key
        for item_compared in list_supplied:
            common_atoms = set(item) & set(item_compared)
            if len(common_atoms) == 2:
                dictionary_results[str(item)] = dictionary_results[str(item)] + 1 #counts how many rings the given ring is fused to
    #print(dictionary_results)

    for value in dictionary_results.values(): #for the dictionary in which the data has been collected, for each ring, we count it to be fused once if at least one common bond has been counted.
        if int(value) >= 1:
            fuses = fuses + 1
    return fuses

#defining functions to test for the following filters - 
#lipinski's, lead-likeness, bioavailablity, ghose, veber, rule of 3, REOS, drug-like QED, and one for testing if all filters were passed or not. for the given file of compounds, none of them pass all the filters.
def lipinski_scorer(mass: float, logP: float, donorCount: int, acceptorCount: int):
    score = 0
    if mass <= 500:
        score = score + 1
    if logP <= 5:
        score = score + 1
    if donorCount <= 5:
        score = score + 1
    if acceptorCount <= 10:
        score = score + 1

    return score

def lead_likeness_scorer(mass: float,logD: float, ringCount: int,rotatableBondCount: int,donorCount: int, acceptorCount: int):
    score = 0
    if mass <= 450:
        score = score + 1
    if logD >= (-4) and logD <= 4:
        score = score + 1
    if ringCount <=  4:
        score = score + 1
    if rotatableBondCount <= 10:
        score = score + 1
    if donorCount <= 5:
        score = score + 1
    if acceptorCount <= 8:
        score = score + 1
    
    return score

def bioavailability_scorer(mass: float, logP: float, donorCount: int, acceptorCount: int, rotatableBondCount: int, PSA: float, fusedAromaticRingCount: int):
    score = 0
    if mass <= 500:
        score = score + 1
    if logP <= 5:
        score = score + 1
    if donorCount <= 5:
        score = score + 1
    if acceptorCount <= 10:
        score = score + 1
    if rotatableBondCount <= 10:
        score = score + 1
    if PSA <= 200:
        score = score + 1
    if fusedAromaticRingCount <= 5:
        score = score + 1

    return score

def ghose_scorer(mass: float, logP: float, atom_count: int, molar_refractivity: float):
    score = 0
    if mass > 160 and mass < 480:
        score = score + 1
    if logP > (-0.4) and logP < 5.6:
        score = score + 1
    if atom_count > 20 and atom_count < 70:
        score = score + 1
    if molar_refractivity > 40 and molar_refractivity < 130:
        score = score + 1

    return score

def veber_scorer(rotatableBondCount: int,TPSA: float):
    score = 0
    if rotatableBondCount <= 10:
        score = score + 1
    if TPSA <= 140:
        score = score + 1

    return score

def REOS_scorer(mass: float, logP: float, donorCount: int, acceptorCount: int, formal_charge: float, rotatableBondCount: int, heavy_atom_count: int):
    score = 0
    if mass > 200 and mass < 500:
        score = score + 1
    if logP > (-5) and logP < 5:
        score = score + 1
    if donorCount > 0 and donorCount < 5:
        score = score + 1
    if acceptorCount > 0 and acceptorCount < 10:
        score = score + 1
    if formal_charge > (-2) and formal_charge < 2:
        score = score + 1
    if rotatableBondCount > 0 and rotatableBondCount < 8:
        score = score + 1
    if heavy_atom_count > 15 and heavy_atom_count < 50:
        score = score + 1

    return score

def rule_of_3_scorer(mass: float, logP: float, donorCount: int, acceptorCount: int, rotatableBondCount: int):
    score = 0
    if mass <= 300:
        score = score + 1
    if logP <= 3:
        score = score + 1
    if donorCount <= 3:
        score = score + 1
    if acceptorCount <= 3:
        score = score + 1
    if rotatableBondCount <= 3:
        score = score + 1

    return score

def drug_like_QED_scorer(mass: float, ringCount: int, rotatableBondCount: int, donorCount: int, acceptorCount: int, logP: float):
    score = 0
    if mass < 400:
        score = score + 1
    if ringCount > 0:
        score = score + 1
    if rotatableBondCount < 5:
        score = score + 1
    if donorCount <= 5:
        score = score + 1
    if acceptorCount <= 10:
        score = score + 1
    if logP < 5:
        score = score + 1

    return score

def all_filters_test(lipinski_score,lead_likeness_score,bioavailability_score,ghose_score,veber_score,REOS_score,rule_of_3_score,drug_like_QED_score):
    if lipinski_score == 4 and lead_likeness_score == 6 and bioavailability_score >= 6 and ghose_score == 4 and veber_score == 2 and REOS_score == 7 and rule_of_3_score == 5 and drug_like_QED_score == 6:
        return True
    else:
        return False

'''
test = 0
for compound in db:
    print(compound.GetNumAtoms()) #Return the number of atoms
    print(Chem.MolToSmiles(compound)) #get the SMILES code
    print(AllChem.Compute2DCoords(compound)) #2D coordinates #all are zero for me
    test = test + 1
print(f"Number of compounds is {test}.") #ok this works

for compound in db:
    print(compound.GetBonds()[0].GetBondType()) #gets bond type of the first bond of the first compound
    for atom in compound.GetAtoms():
        print(atom.GetAtomicNum()) #gets the atomic number
        break
'''
'''        
for compound in db:
    tmp=AllChem.Compute2DCoords(compound)
'''

#this loop is for my personal use. it extract the dictionary of the descriptors of the first compound.
for compound in db:
    descriptors_list = []
    for key in (Descriptors.CalcMolDescriptors(compound)).keys():
        descriptors_list.append(key)
    break

#creating the database schema. first creating a connection and cursor.
connection = sqlite3.connect(database_path)
cur = connection.cursor()

#making sure the ddl file exists
try:
    assert os.path.isfile('sqlddl.sql') == True
except AssertionError:
    logger.error(f"Hello, the sql ddl script was not found in the current working directory. Please fix it and re-run the script. Raising SystemExit.")
    raise SystemExit(1)

#reading the sql script.
with open('sqlddl.sql') as sqlfile:
    sql_script = sqlfile.read() #reading the script and saving it

#print(sql_script)
cur.executescript(sql_script) #executing the script in the given schema.
connection.commit() #committing the changes to the database
cur.close()
connection.close() #closing the cursor and the connection.

#opening another connection to load the database
connection = sqlite3.connect(database_path)
cur = connection.cursor()

#this extracts the parameters for the database.
for compound in db:
    descriptors = Descriptors.CalcMolDescriptors(compound) #this returns a dictionary, where the descriptor is the key and the value is the respective value.

    #QED stands for quantitative estimation of drug-likeness
    qed = descriptors['qed']
    #print(qed)
    qed = round(qed,3)
    #print(qed)

    #We propose the spacial score (SPS) as an empirical scoring system that builds upon the principle underlying Fsp3 and FCstereo and expresses the spacial complexity of a compound in a uniform manner on a highly granular scale.
    SPS = descriptors['SPS']
    #print(SPS)
    SPS = round(SPS,3)
    #print(SPS)
    
    '''
    Exact Molecular Mass versus Molecular Weight
    "In the discussion in the last section, the term exact mass was introduced. Again this is the mass of a molecule calculated with only the most abundant isotopes present; these are usually the lightest isotopes. The molecular weights given in the Handbook of Chemistry and Physics are different; they are the average molecular weights for all the molecules with different isotopic compositions present in a compound, each weighted for the natural abundance of the respective isotope or isotopes. For example the molecular weight of propane is 44.0956 g/mol and for acetaldehyde, 44.0526 g/mol. The exact molecular masses for the molecules with the lightest isotopes are 44.0626 amu and 44.0262 amu, respectively. Note, we use different units for molecular weight than we use for molecular mass."
    Copyright information: Original content © University of Colorado, Boulder, Chemistry and Biochemistry Department, 2011. The information on these pages is available for academic use without restriction.
    Accessed: https://www.orgchemboulder.com/Spectroscopy/MS/molmassmw.shtml [15 February 2024]
    '''

    MolWt = descriptors['MolWt']
    MolWt = round(MolWt,3)
    #print(MolWt)
    ExactMolWt = descriptors['ExactMolWt']
    ExactMolWt = round(ExactMolWt,3)
    #print(ExactMolWt)

    NumValenceElectrons = descriptors['NumValenceElectrons']
    #print(NumValenceElectrons)

    NumRadicalElectrons = descriptors['NumRadicalElectrons']
    #print(NumRadicalElectrons)

    HBA = descriptors['NOCount']
    #print(HBA)
    #print(f"flag for me {descriptors['NumHAcceptors']}")

    HBD = descriptors['NHOHCount']
    #print(HBD)
    #print(f"flag for me {descriptors['NumHDonors']}")
    #those two flags returned the same values

    LogP = descriptors['MolLogP']
    LogP = round(LogP,3)
    #print(LogP)

    NumAromaticRings = descriptors['NumAromaticRings']
    #print(NumAromaticRings)

    #print(descriptors)
    #print('flag 1')

    #getting the formula
    formula = Chem.rdMolDescriptors.CalcMolFormula(compound)
    #print(formula)

    name = compound.GetProp('Name')
    #print(name)

    LogD = compound.GetProp('LogD')
    LogD = round(float(LogD),3)
    #print(LogD)

    RingCount = descriptors['RingCount']
    #print(RingCount)
    #print('flag 2')

    NumRotatableBonds = descriptors['NumRotatableBonds']
    #print(NumRotatableBonds)

    #something abour surface area
    TPSA = descriptors['TPSA']
    TPSA = round(TPSA,3)
    #print(TPSA)

    molar_refractivity = Chem.Crippen.MolMR(compound)
    molar_refractivity = round(molar_refractivity,3)
    #print(molar_refractivity)

    Mol_ID = compound.GetProp('Mol_ID')
    #print(Mol_ID)

    CdId = compound.GetProp('CdId')
    #print(CdId)

    number_of_atoms = Chem.rdchem.Mol.GetNumAtoms(compound)
    #print(number_of_atoms)

    formal_charge = Chem.rdmolops.GetFormalCharge(compound)
    #print(formal_charge)

    heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(compound)
    #print(heavy_atoms) #an isotopic form of an atom that contains more than the common number of neutrons and is thus of greater relative atomic mass than the most abundant or most commonly observed isotope.

    #this extracts a list of lists where each list contains information about a ring that is aromatic. the said information is a tuple of atom numbers contained in that ring (with respect to the whole compound, so numbering/labelling is consistent.)
    #for example, this is the list for the first compound (CdId 1):
    #[(4, 5, 6, 15, 16, 3), (7, 8, 9, 14, 5, 6), (10, 11, 12, 13, 14, 9)]
    #this compound has three aromatic rings, so the list contains three tuples. each tuple says which atoms are in a ring. therefore, the first compound's atom numbers (4, 5, 6, 15, 16, 3) are in an aromatic ring.
    aromatic_rings = [ring for ring in compound.GetRingInfo().AtomRings() if all(compound.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring)]
    #print(aromatic_rings)

    fused_aro_rings = fused_aro_rings_calculator(aromatic_rings)
    #print(fused_aro_rings)

    lipinski_score = lipinski_scorer(mass=float(MolWt),logP=float(LogP),donorCount=int(HBD),acceptorCount=int(HBA))
    #print(lipinski_score)
    #print('flag 3')

    lead_likeness_score = lead_likeness_scorer(mass=float(MolWt),logD=float(LogD),ringCount=int(RingCount),rotatableBondCount=int(NumRotatableBonds),donorCount=int(HBD),acceptorCount=int(HBA))
    #print(lead_likeness_score)

    bioavailability_score = bioavailability_scorer(mass=float(MolWt),logP=float(LogP),donorCount=int(HBD),acceptorCount=int(HBA),rotatableBondCount=int(NumRotatableBonds),PSA=float(TPSA),fusedAromaticRingCount=int(fused_aro_rings))
    #print(bioavailability_score)

    ghose_score = ghose_scorer(float(MolWt),float(LogP),int(number_of_atoms),float(molar_refractivity))
    #print(ghose_score)

    veber_score = veber_scorer(int(NumRotatableBonds),float(TPSA))
    #print(veber_score)

    REOS_score = REOS_scorer(float(MolWt),float(LogP),int(HBD),int(HBA),float(formal_charge),int(NumRotatableBonds),int(heavy_atoms))
    #print(REOS_score)

    rule_of_3_score = rule_of_3_scorer(float(MolWt),float(LogP),int(HBD),int(HBA),int(NumRotatableBonds))
    #print(rule_of_3_score)

    drug_like_QED_score = drug_like_QED_scorer(float(MolWt),int(RingCount),int(NumRotatableBonds),int(HBD),int(HBA),float(LogP))
    #print(drug_like_QED_score)

    smiles = Chem.MolToSmiles(compound)
    #print(smiles)

    try:
        IUPAC = translate_forward(smiles)
    except: #translate_forward is a function from the library STOUT which is a neural net trained on many known SMILES and their respective IUPAC codes.
        IUPAC = None
    #print(IUPAC)
    #print(name)

    #print('flag 4')

    Draw.MolToFile(db[(int(CdId)-1)],'image.png') #drawing the image using the CdId as an index to access the original database created from the sdf file.

    with open('image.png', 'rb') as image_file:
        binaryData = image_file.read() #converting the image into binary data (BLOB) so it can be inserted into the data using SQL.
        
    #print(binaryData)

    all_filters_boolean = all_filters_test(lipinski_score,lead_likeness_score,bioavailability_score,ghose_score,veber_score,REOS_score,rule_of_3_score,drug_like_QED_score)

    #now the instance of data to be input into the database can be created in the form of a list of parameters. #this follows the chronological order in which the columns exist in the sql ddl file, so it respects the database schema.
    list_to_input = [CdId,name,formula,binaryData,IUPAC,Mol_ID,MolWt,ExactMolWt,LogP,LogD,HBD,HBA,NumRotatableBonds,RingCount,TPSA,NumAromaticRings,fused_aro_rings,smiles,qed,SPS,NumValenceElectrons,NumRadicalElectrons,molar_refractivity,number_of_atoms,formal_charge,heavy_atoms,lipinski_score,lead_likeness_score,bioavailability_score,ghose_score,veber_score,REOS_score,rule_of_3_score,drug_like_QED_score,all_filters_boolean]
    #print(list_to_input)
    #print(len(list_to_input)) #we have 35 attributes of the entity Compounds.

    #creating the insert statement
    question_marks = '?,' * 35
    question_marks = question_marks.rstrip(',')
    insert_statement = "INSERT INTO Compounds VALUES (%s);" % question_marks
    #print(question_marks)
    #print(insert_statement)

    #the data is now ready to be inserted.
    try:
        cur.execute(insert_statement,list_to_input)
        connection.commit() #committing the changes to the database
    except:
        logger.warning(f"The data for compound {CdId} could not be inserted into the database, skipping this.")
        continue

#print(descriptors_list)

#closing the filehandle to the database
cur.close()
connection.close()

#removing the image file
os.remove("image.png")

#for 100 compounds, this takes about 1.5 hours to run on a regular virtual machine.