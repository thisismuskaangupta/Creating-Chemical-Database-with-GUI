'''
AI-tools like GPT and internet forums like StackOverflow were used to produce parts of this script, particularly harder functionality like sorting columns. Code was adapted and customized to fulfil the needs of this assignment. 
'''

#importing logging
try:
    import logging
except ModuleNotFoundError:
    print("Hello, the logging module is not installed. Please install it and run the script again. Raising SystemExit.")
    raise SystemExit(1)

#setting up the logger
#getting the logger
loggerA = logging.getLogger('logger for GUI script')
#setting the level of the logger
loggerA.setLevel(logging.WARNING)
#creating a stream handler.
shandler = logging.StreamHandler()
#formatting the stream handler
shandler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - {%(pathname)s:%(lineno)d} - %(message)s'))
#adding the stream handler to the logger.
loggerA.addHandler(shandler)

#importing packages
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import Descriptors
    import sqlite3
    import tkinter as tk
    from tkinter import ttk
    from tkinter import PhotoImage
    from PIL import Image, ImageTk
    import io
    import os
except ModuleNotFoundError as e:
    loggerA.error(f"Hello, the following module was not found, please download it and run the script again. Raising systemexit.\n{e}")
    raise SystemExit(1)

#setting file paths
database_path = '2873826.db'
sdf_filepath = 'Molecules8.sdf'

#making sure the sdf database filepath exists
try:
    assert os.path.isfile(sdf_filepath) == True
except AssertionError:
    loggerA.error(f"Hello the given sdf file {sdf_filepath} does not exist. Please fix the error or make sure the file exists. Raising SystemExit.")
    raise SystemExit(1)

#this function is adapted from AI-tools and several StackOverflow posts.
#defining a function to sort the columns in the GUI by clicking on them.
def treeview_sort_column(tree, col, reverse):
    # extracting all the data from the columns of the TreeView tree and storing it in a list with the item ID and value
    data = [(tree.set(k, col), k) for k in tree.get_children('')]
    
    # converting values to floats for numeric comparison
    try:
        data.sort(key=lambda t: float(t[0]), reverse=reverse)
    except ValueError: #if they are not numbers and so cannot be converted to floats, string comparison is done
        data.sort(key=str, reverse=reverse)
    
    # rearranging the items in the sorted order
    for index, (_, k) in enumerate(data):
        tree.move(k, '', index)

    # implementing functionality to make sure that if the column name heading is clicked again, sorting is done in reverse order. the default is increasing order.
    tree.heading(col, command=lambda _col=col: treeview_sort_column(tree, _col, not reverse))

#defining functions to be called whenever buttons in the GUI are clicked.
def button1clicker(): #this is lipinski
    global tree, window
    for i in tree.get_children():
        #print(i)
        #print((tree.item(i)))
        if ((tree.item(i)['values'])[26]) != 4: #the lipinski score is checked, if it is not equal to 4 (which is when the compound passes the lipinski test), the row is deleted.
            tree.delete(i)
    window.update() #the window is updated to reflect the changes.

#all the filters work by a similar logic
def button2clicker(): #this is veber
    global tree, window
    for i in tree.get_children():
        #print(i)
        #print((tree.item(i)))
        if ((tree.item(i)['values'])[30]) != 2:
            tree.delete(i)
    window.update()

def button3clicker(): #this is lead likeness
    global tree, window
    for i in tree.get_children():
        #print(i)
        #print((tree.item(i)))
        if ((tree.item(i)['values'])[27]) != 6:
            tree.delete(i)
    window.update()

def button4clicker(): #this is REOS
    global tree, window
    for i in tree.get_children():
        #print(i)
        #print((tree.item(i)))
        if ((tree.item(i)['values'])[31]) != 7:
            tree.delete(i)
    window.update()

def button5clicker(): #this is bioavailability
    global tree, window
    for i in tree.get_children():
        #print(i)
        #print((tree.item(i)))
        if ((tree.item(i)['values'])[28]) < 6:
            tree.delete(i)
    window.update()

def button6clicker(): #this is rule of 3
    global tree, window
    for i in tree.get_children():
        #print(i)
        #print((tree.item(i)))
        if ((tree.item(i)['values'])[32]) != 5:
            tree.delete(i)
    window.update()

def button7clicker(): #this is ghose
    global tree, window
    for i in tree.get_children():
        #print(i)
        #print((tree.item(i)))
        if ((tree.item(i)['values'])[29]) != 4:
            tree.delete(i)
    window.update()

def button8clicker(): #this is drug like qed
    global tree, window
    for i in tree.get_children():
        #print(i)
        #print((tree.item(i)))
        if ((tree.item(i)['values'])[33]) != 6:
            tree.delete(i)
    window.update()

def button9clicker(): #this is the clear button, upon clicking this button, we want the entry boxes to be cleared and the original database to be loaded up.

    global cur1, tree, window

    for i in tree.get_children():
        tree.delete(i)
    window.update() #first the entire tree is deleted so the original tree can be created again.

    cur1.execute("SELECT * FROM Compounds") #the entire database is selected and extracted
    rows = cur1.fetchall()  


    if not hasattr(tree, 'image_cache'): #an empty list for caching the images is created
        tree.image_cache = []

    for row in rows:
        # Getting the BLOB data for the image from the fourth column (according to 1 indexing)
        image_data = row[3]
        
        # Converting the BLOB data to an image
        try:
            image_extracted = Image.open(io.BytesIO(image_data)) #opening image
            image_extracted.thumbnail((200,200))  # Resizing the image to fit the GUI
            photo_to_insert = ImageTk.PhotoImage(image_extracted)  
        except Exception as e:
            loggerA.warning(f"Hello, the image for the compound {row[0]} could not be loaded, skipping this line.")
            continue  
        
        # creating the values to insert in the TreeView
        values = row
        
        # Inserting the row into the treeview
        tree.insert("", tk.END, values=values, image=photo_to_insert,text=row[0])
        
        # Caching the object photo_to_insert to prevent it from being garbage collected
        tree.image_cache.append(photo_to_insert)

    #further, deleting all the entries the user may have made in the GUI
    entry1.delete(0,tk.END)
    entry2.delete(0,tk.END)
    entry3.delete(0,tk.END)
    entry4.delete(0,tk.END)
    entry5.delete(0,tk.END)
    entry6.delete(0,tk.END)
    entry7.delete(0,tk.END)
    entry8.delete(0,tk.END)
    entry9.delete(0,tk.END)
    entry10.delete(0,tk.END)
    entry11.delete(0,tk.END)
    entry12.delete(0,tk.END)
    entry13.delete(0,tk.END)
    entry14.delete(0,tk.END)
    entry15.delete(0,tk.END)
    entry16.delete(0,tk.END)
    entry17.delete(0,tk.END)
    entry18.delete(0,tk.END)
    entry19.delete(0,tk.END)
    entry20.delete(0,tk.END)
    entry21.delete(0,tk.END)
    entry22.delete(0,tk.END)

def button10clicker(): #this is the submit button

    #we extract whatever the user has passed in the GUI and use it to overwrite the default
    global min_mol_weight, max_mol_weight, min_logp, max_logp, min_logd, max_logd, min_hbd, max_hbd, min_hba, max_hba, min_rot_bonds, max_rot_bonds, min_ring_count, max_ring_count, min_TPSA, max_TPSA, min_qed, max_qed, min_mr, max_mr, min_no_of_atoms, max_no_of_atoms,cur1,tree,window
    
    if entry1.get() != '':
        min_mol_weight = float(entry1.get())

    if entry2.get() != '':
        min_logp = float(entry2.get())

    if entry3.get() != '':
        min_logd = float(entry3.get())

    if entry4.get() != '':
        min_hbd = float(entry4.get())

    if entry5.get() != '':
        min_hba = float(entry5.get())

    if entry6.get() != '':
        min_rot_bonds = float(entry6.get())

    if entry7.get() != '':
        min_ring_count = float(entry7.get())

    if entry8.get() != '':
        min_TPSA = float(entry8.get())

    if entry9.get() != '':
        min_qed = float(entry9.get())

    if entry10.get() != '':
        min_mr = float(entry10.get())

    if entry11.get() != '':
        min_no_of_atoms = float(entry11.get())

    if entry12.get() != '':
        max_mol_weight = float(entry12.get())

    if entry13.get() != '':
        max_logp = float(entry13.get())

    if entry14.get() != '':
        max_logd = float(entry14.get())

    if entry15.get() != '':
        max_hbd = float(entry15.get())

    if entry16.get() != '':
        max_hba = float(entry16.get())

    if entry17.get() != '':
        max_rot_bonds = float(entry17.get())

    if entry18.get() != '':
        max_ring_count = float(entry18.get())

    if entry19.get() != '':
        max_TPSA = float(entry19.get())

    if entry20.get() != '':
        max_qed = float(entry20.get())

    if entry21.get() != '':
        max_mr = float(entry21.get())

    if entry22.get() != '':
        max_no_of_atoms = float(entry22.get())

    for i in tree.get_children():
        tree.delete(i)
    window.update() #all the entries in the TreeView are deleted so that the selected ones only may be displayed.

    select_statement = "SELECT * FROM Compounds WHERE Molecular_Weight BETWEEN %s AND %s AND LogP BETWEEN %s AND %s AND LogD BETWEEN %s AND %s AND Hydrogen_Bond_Donors BETWEEN %s AND %s AND Hydrogen_Bond_Acceptors BETWEEN %s AND %s AND Rotatable_Bonds BETWEEN %s AND %s AND Number_of_Rings BETWEEN %s AND %s AND TPSA BETWEEN %s AND %s AND QED BETWEEN %s AND %s AND Molar_Refractivity BETWEEN %s AND %s AND Number_of_Atoms BETWEEN %s AND %s;" %(min_mol_weight,max_mol_weight,min_logp,max_logp,min_logd,max_logd,min_hbd,max_hbd,min_hba,max_hba,min_rot_bonds,max_rot_bonds,min_ring_count,max_ring_count,min_TPSA,max_TPSA,min_qed,max_qed,min_mr,max_mr,min_no_of_atoms,max_no_of_atoms) #we use defaults AND user-passed inputs to query the database using a single cumulative select statement.
    rows = cur1.execute(select_statement).fetchall() #the query results are extracted

    #a similar logic as the clear button is used to finally write out the results to the GUI
    if not hasattr(tree, 'image_cache'):
        tree.image_cache = []

    for row in rows:
        image_data = row[3]
        
        # Convert the BLOB data to an image
        try:
            image_extracted = Image.open(io.BytesIO(image_data))
            image_extracted.thumbnail((200,200))  
            photo_to_insert = ImageTk.PhotoImage(image_extracted)  
        except Exception as e:
            loggerA.warning(f"Hello, the image for the compound {row[0]} could not be loaded, skipping this line.")
            continue  
        
        values = row
        
        tree.insert("", tk.END, values=values, image=photo_to_insert,text=row[0])
        
        tree.image_cache.append(photo_to_insert)
    
    window.update()

#defining a function to export results into a csv file
def button11clicker():
    global tree
    with open('results.csv','w') as out_file:
        for i in tree.get_children():
            out_file.write(f"{str(tree.item(i)['values'][:2])}\n")

#creating the gui

#creating a window
window = tk.Tk()
greeting = tk.Label(text="Chemical Database") #writing a heading for the whole GUI
greeting.pack()

#creating two frames
frame2 = tk.Frame(master=window, width=400, height=1000, bg="PaleVioletRed1")
frame2.pack(fill=tk.BOTH, side=tk.LEFT,expand=True)

frame3 = tk.Frame(master=window, width=800, bg="yellow")
frame3.pack(fill=tk.BOTH, side=tk.RIGHT,expand=True)

#creating column names, we have 35 columns and so we can use a loop
column_names = []
for i in range(1,36,):
    name = 'c' + str(i)
    column_names.append(name)

s=ttk.Style() #using a style to customize the row height so that the images may be added conveniently
s.theme_use('clam')
s.configure('Treeview', rowheight=200)

#creating the tree
tree = ttk.Treeview(frame3, column=column_names, show='tree headings')

#adding the functionality to sort the column by clicking on it
for col in column_names:
    tree.heading(col, text=col, command=lambda _col=col: treeview_sort_column(tree, _col, False))

#naming all the columns in the GUI tree
tree.column("#1", anchor=tk.CENTER,stretch=False)
tree.heading("#1", text="Candidate_ID")
tree.column("#2", anchor=tk.CENTER,stretch=False)
tree.heading("#2", text="Molecule_Name")
tree.column("#3", anchor=tk.CENTER,stretch=False)
tree.heading("#3", text="Formula")
tree.column("#4", anchor=tk.CENTER,stretch=False,width=0,minwidth=0) #hiding the images column which will display BLOB data, which is not readable by the human eye and the images are displayed in a different 0th column.
tree.heading("#4", text="Structure")
tree.column("#5", anchor=tk.CENTER,stretch=False)
tree.heading("#5", text="IUPAC")
tree.column("#6", anchor=tk.CENTER,stretch=False)
tree.heading("#6", text="Molecule_Number")
tree.column("#7", anchor=tk.CENTER,stretch=False)
tree.heading("#7", text="Molecular_Weight")
tree.column("#8", anchor=tk.CENTER,stretch=False)
tree.heading("#8", text="Exact_Molecular_Weight")
tree.column("#9", anchor=tk.CENTER,stretch=False)
tree.heading("#9", text="LogP")
tree.column("#10", anchor=tk.CENTER,stretch=False)
tree.heading("#10", text="LogD")
tree.column("#11", anchor=tk.CENTER,stretch=False)
tree.heading("#11", text="Hydrogen_Bond_Donors")
tree.column("#12", anchor=tk.CENTER,stretch=False)
tree.heading("#12", text="Hydrogen_Bond_Acceptors")
tree.column("#13", anchor=tk.CENTER,stretch=False)
tree.heading("#13", text="Rotatable_Bonds")
tree.column("#14", anchor=tk.CENTER,stretch=False)
tree.heading("#14", text="Number_of_Rings")
tree.column("#15", anchor=tk.CENTER,stretch=False)
tree.heading("#15", text="TPSA")
tree.column("#16", anchor=tk.CENTER,stretch=False)
tree.heading("#16", text="Number_of_Aromatic_Rings")
tree.column("#17", anchor=tk.CENTER,stretch=False)
tree.heading("#17", text="Number_of_Fused_Aromatic_Rings")
tree.column("#18", anchor=tk.CENTER,stretch=False)
tree.heading("#18", text="SMILES")
tree.column("#19", anchor=tk.CENTER,stretch=False)
tree.heading("#19", text="QED")
tree.column("#20", anchor=tk.CENTER,stretch=False)
tree.heading("#20", text="SPS")
tree.column("#21", anchor=tk.CENTER,stretch=False)
tree.heading("#21", text="Valence_Electrons")
tree.column("#22", anchor=tk.CENTER,stretch=False)
tree.heading("#22", text="Radical_Electrons")
tree.column("#23", anchor=tk.CENTER,stretch=False)
tree.heading("#23", text="Molar_Refractivity")
tree.column("#24", anchor=tk.CENTER,stretch=False)
tree.heading("#24", text="Number_of_Atoms")
tree.column("#25", anchor=tk.CENTER,stretch=False)
tree.heading("#25", text="Formal_Charge")
tree.column("#26", anchor=tk.CENTER,stretch=False)
tree.heading("#26", text="Number_of_Heavy_Atoms")
tree.column("#27", anchor=tk.CENTER,stretch=False)
tree.heading("#27", text="Lipinski_Score")
tree.column("#28", anchor=tk.CENTER,stretch=False)
tree.heading("#28", text="Lead_Likeness_Score")
tree.column("#29", anchor=tk.CENTER,stretch=False)
tree.heading("#29", text="Bioavailability_Score")
tree.column("#30", anchor=tk.CENTER,stretch=False)
tree.heading("#30", text="Ghose_Score")
tree.column("#31", anchor=tk.CENTER,stretch=False)
tree.heading("#31", text="Veber_Score")
tree.column("#32", anchor=tk.CENTER,stretch=False)
tree.heading("#32", text="REOS_Score")
tree.column("#33", anchor=tk.CENTER,stretch=False)
tree.heading("#33", text="Rule_of_3_Score")
tree.column("#34", anchor=tk.CENTER,stretch=False)
tree.heading("#34", text="Drug_Like_QED_Score")
tree.column("#35", anchor=tk.CENTER,stretch=False)
tree.heading("#35", text="All_Filters_Test")

#adding scroll bars
scrollbar_x = ttk.Scrollbar(frame3, command=tree.xview, orient=tk.HORIZONTAL)
scrollbar_y = ttk.Scrollbar(frame3,command=tree.yview,orient=tk.VERTICAL)
tree.config(xscrollcommand=scrollbar_x.set,yscrollcommand=scrollbar_y.set)
scrollbar_x.pack(fill=tk.X, side=tk.BOTTOM)
scrollbar_y.pack(fill=tk.Y,side=tk.RIGHT)

tree.pack(fill=tk.BOTH, expand=True) #packing all the code into the tree

#creating a connection and a cursor to the database
con1 = sqlite3.connect(database_path)
cur1 = con1.cursor()
#the default view of the tree is loaded, the default view being the entire database
cur1.execute("SELECT * FROM Compounds")
rows = cur1.fetchall()  

#just as before, the tree is populated
if not hasattr(tree, 'image_cache'):
    tree.image_cache = []

for row in rows:
    image_data = row[3]
    
    try:
        image_extracted = Image.open(io.BytesIO(image_data))
        image_extracted.thumbnail((200,200))  
        photo_to_insert = ImageTk.PhotoImage(image_extracted)  
    except Exception as e:
        loggerA.warning(f"Hello, the image for the compound {row[0]} could not be loaded, skipping this line.")
        continue  
        
    values = row
    
    tree.insert("", tk.END, values=values, image=photo_to_insert,text=row[0])
    
    tree.image_cache.append(photo_to_insert)

tree['displaycolumns'] = tuple(range(len(column_names))) #making sure the display is executed correctly

#setting default settings, by writing SQL to extract minimum and maximum of all the chosen parameters and subtracting or adding one to get the whole range. these defaults are used to capture the entire database, unless the user chooses to filter it.
min_mol_weight = cur1.execute("SELECT MIN(Molecular_Weight) FROM Compounds;").fetchall()
min_mol_weight = float((min_mol_weight[0])[0]) - 1
#print(min_mol_weight)

max_mol_weight = cur1.execute("SELECT MAX(Molecular_Weight) FROM Compounds;").fetchall()
max_mol_weight = float((max_mol_weight[0])[0]) + 1
#print(max_mol_weight)

min_logp = cur1.execute("SELECT MIN(LogP) FROM Compounds;").fetchall()
min_logp = float((min_logp[0])[0]) - 1

max_logp = cur1.execute("SELECT MAX(LogP) FROM Compounds;").fetchall()
max_logp = float((max_logp[0])[0]) + 1

min_logd = cur1.execute("SELECT MIN(LogD) FROM Compounds;").fetchall()
min_logd = float((min_logd[0])[0]) - 1

max_logd = cur1.execute("SELECT MAX(LogD) FROM Compounds;").fetchall()
max_logd = float((max_logd[0])[0]) + 1

min_hbd = cur1.execute("SELECT MIN(Hydrogen_Bond_Donors) FROM Compounds;").fetchall()
min_hbd = float((min_hbd[0])[0]) - 1

max_hbd = cur1.execute("SELECT MAX(Hydrogen_Bond_Donors) FROM Compounds;").fetchall()
max_hbd = float((max_hbd[0])[0]) + 1

min_hba = cur1.execute("SELECT MIN(Hydrogen_Bond_Acceptors) FROM Compounds;").fetchall()
min_hba = float((min_hba[0])[0]) - 1

max_hba = cur1.execute("SELECT MAX(Hydrogen_Bond_Acceptors) FROM Compounds;").fetchall()
max_hba = float((max_hba[0])[0]) + 1

min_rot_bonds = cur1.execute("SELECT MIN(Rotatable_Bonds) FROM Compounds;").fetchall()
min_rot_bonds = float((min_rot_bonds[0])[0]) - 1

max_rot_bonds = cur1.execute("SELECT MAX(Rotatable_Bonds) FROM Compounds;").fetchall()
max_rot_bonds = float((max_rot_bonds[0])[0]) + 1

min_ring_count = cur1.execute("SELECT MIN(Number_of_Rings) FROM Compounds;").fetchall()
min_ring_count = float((min_ring_count[0])[0]) - 1

max_ring_count = cur1.execute("SELECT MAX(Number_of_Rings) FROM Compounds;").fetchall()
max_ring_count = float((max_ring_count[0])[0]) + 1

min_TPSA = cur1.execute("SELECT MIN(TPSA) FROM Compounds;").fetchall()
min_TPSA = float((min_TPSA[0])[0]) - 1

max_TPSA = cur1.execute("SELECT MAX(TPSA) FROM Compounds;").fetchall()
max_TPSA = float((max_TPSA[0])[0]) + 1

min_qed = cur1.execute("SELECT MIN(QED) FROM Compounds;").fetchall()
min_qed = float((min_qed[0])[0]) - 1

max_qed = cur1.execute("SELECT MAX(QED) FROM Compounds;").fetchall()
max_qed = float((max_qed[0])[0]) + 1

min_mr = cur1.execute("SELECT MIN(Molar_Refractivity) FROM Compounds;").fetchall()
min_mr = float((min_mr[0])[0]) - 1

max_mr = cur1.execute("SELECT MAX(Molar_Refractivity) FROM Compounds;").fetchall()
max_mr = float((max_mr[0])[0]) + 1

min_no_of_atoms = cur1.execute("SELECT MIN(Number_of_Atoms) FROM Compounds;").fetchall()
min_no_of_atoms = float((min_no_of_atoms[0])[0]) - 1

max_no_of_atoms = cur1.execute("SELECT MAX(Number_of_Atoms) FROM Compounds;").fetchall()
max_no_of_atoms = float((max_no_of_atoms[0])[0]) + 1

#adding labels
label1 = tk.Label(master=frame2, text="Parameters")
label1.place(relx=0.5, rely=0.025, anchor=tk.CENTER)

label2 = tk.Label(master=frame2,text="From")
label2.place(relx=0.45,rely=0.075,anchor='center')

label3 = tk.Label(master=frame2,text='To')
label3.place(relx=0.75,rely=0.075,anchor='center')

label4 = tk.Label(master=frame2,text='Molecular Weight')
label4.place(relx=0.15,rely=0.125,anchor='center')

#alongwith the labels, adding text boxes to collect data the user wishes to input
entry1 = tk.Entry(bg="white",master=frame2)
entry1.place(anchor="center",relx=0.45,rely=0.125,width=100)
entry12 = tk.Entry(bg="white",master=frame2)
entry12.place(anchor="center",relx=0.75,rely=0.125,width=100)

label5 = tk.Label(master=frame2,text='LogP')
label5.place(relx=0.15,rely=0.165,anchor='center')

entry2 = tk.Entry(bg="white",master=frame2)
entry2.place(anchor="center",relx=0.45,rely=0.165,width=100)
entry13 = tk.Entry(bg="white",master=frame2)
entry13.place(anchor="center",relx=0.75,rely=0.165,width=100)

label6 = tk.Label(master=frame2,text='LogD')
label6.place(relx=0.15,rely=0.205,anchor='center')

entry3 = tk.Entry(bg="white",master=frame2)
entry3.place(anchor="center",relx=0.45,rely=0.205,width=100)
entry14 = tk.Entry(bg="white",master=frame2)
entry14.place(anchor="center",relx=0.75,rely=0.205,width=100)

label7 = tk.Label(master=frame2,text='H Bond Donors')
label7.place(relx=0.15,rely=0.245,anchor='center')

entry4 = tk.Entry(bg="white",master=frame2)
entry4.place(anchor="center",relx=0.45,rely=0.245,width=100)
entry15 = tk.Entry(bg="white",master=frame2)
entry15.place(anchor="center",relx=0.75,rely=0.245,width=100)

label8 = tk.Label(master=frame2,text='H Bond Acceptors')
label8.place(relx=0.15,rely=0.285,anchor='center')

entry5 = tk.Entry(bg="white",master=frame2)
entry5.place(anchor="center",relx=0.45,rely=0.285,width=100)
entry16 = tk.Entry(bg="white",master=frame2)
entry16.place(anchor="center",relx=0.75,rely=0.285,width=100)

label9 = tk.Label(master=frame2,text='Rotatable Bonds')
label9.place(relx=0.15,rely=0.325,anchor='center')

entry6 = tk.Entry(bg="white",master=frame2)
entry6.place(anchor="center",relx=0.45,rely=0.325,width=100)
entry17 = tk.Entry(bg="white",master=frame2)
entry17.place(anchor="center",relx=0.75,rely=0.325,width=100)

label10 = tk.Label(master=frame2,text='Ring Count')
label10.place(relx=0.15,rely=0.365,anchor='center')

entry7 = tk.Entry(bg="white",master=frame2)
entry7.place(anchor="center",relx=0.45,rely=0.365,width=100)
entry18 = tk.Entry(bg="white",master=frame2)
entry18.place(anchor="center",relx=0.75,rely=0.365,width=100)

label11 = tk.Label(master=frame2,text='TPSA')
label11.place(relx=0.15,rely=0.405,anchor='center')

entry8 = tk.Entry(bg="white",master=frame2)
entry8.place(anchor="center",relx=0.45,rely=0.405,width=100)
entry19 = tk.Entry(bg="white",master=frame2)
entry19.place(anchor="center",relx=0.75,rely=0.405,width=100)

label12 = tk.Label(master=frame2,text='QED')
label12.place(relx=0.15,rely=0.445,anchor='center')

entry9 = tk.Entry(bg="white",master=frame2)
entry9.place(anchor="center",relx=0.45,rely=0.445,width=100)
entry20 = tk.Entry(bg="white",master=frame2)
entry20.place(anchor="center",relx=0.75,rely=0.445,width=100)

label13 = tk.Label(master=frame2,text='Molar Refractivity')
label13.place(relx=0.15,rely=0.485,anchor='center')

entry10 = tk.Entry(bg="white",master=frame2)
entry10.place(anchor="center",relx=0.45,rely=0.485,width=100)
entry21 = tk.Entry(bg="white",master=frame2)
entry21.place(anchor="center",relx=0.75,rely=0.485,width=100)

label14 = tk.Label(master=frame2,text='No. of Atoms')
label14.place(relx=0.15,rely=0.525,anchor='center')

entry11 = tk.Entry(bg="white",master=frame2)
entry11.place(anchor="center",relx=0.45,rely=0.525,width=100)
entry22 = tk.Entry(bg="white",master=frame2)
entry22.place(anchor="center",relx=0.75,rely=0.525,width=100)

label15 = tk.Label(master=frame2, text="Filters")
label15.place(relx=0.5, rely=0.6, anchor=tk.CENTER)

#adding all the buttons for the filters, and the submit and clear buttons
button1 = tk.Button(frame2, text = 'Lipinski',
                      bg = 'DeepPink4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DeepPink4", 
                      highlightbackground="DeepPink4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button1clicker)
button1.place(relx=0.25,rely=0.6575,anchor=tk.CENTER)

button2 = tk.Button(frame2, text = 'Veber',
                      bg = 'DeepPink4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DeepPink4", 
                      highlightbackground="DeepPink4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button2clicker)
button2.place(relx=0.75,rely=0.6575,anchor=tk.CENTER)

button3 = tk.Button(frame2, text = 'Lead Likeness',
                      bg = 'DeepPink4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DeepPink4", 
                      highlightbackground="DeepPink4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button3clicker)
button3.place(relx=0.25,rely=0.7325,anchor=tk.CENTER)

button4 = tk.Button(frame2, text = 'REOS',
                      bg = 'DeepPink4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DeepPink4", 
                      highlightbackground="DeepPink4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button4clicker)
button4.place(relx=0.75,rely=0.7325,anchor=tk.CENTER)

button5 = tk.Button(frame2, text = 'Bioavailability',
                      bg = 'DeepPink4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DeepPink4", 
                      highlightbackground="DeepPink4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button5clicker)
button5.place(relx=0.25,rely=0.8075,anchor=tk.CENTER)

button6 = tk.Button(frame2, text = 'Rule of 3',
                      bg = 'DeepPink4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DeepPink4", 
                      highlightbackground="DeepPink4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button6clicker)
button6.place(relx=0.75,rely=0.8075,anchor=tk.CENTER)

button7 = tk.Button(frame2, text = 'Ghose',
                      bg = 'DeepPink4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DeepPink4", 
                      highlightbackground="DeepPink4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button7clicker)
button7.place(relx=0.25,rely=0.8825,anchor=tk.CENTER)

button8 = tk.Button(frame2, text = 'Drug-Like QED',
                      bg = 'DeepPink4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DeepPink4", 
                      highlightbackground="DeepPink4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button8clicker)
button8.place(relx=0.75,rely=0.8825,anchor=tk.CENTER)

button9 = tk.Button(frame2, text = 'Clear',
                      bg = 'DarkOrchid4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DarkOrchid4", 
                      highlightbackground="DarkOrchid4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button9clicker)
button9.place(relx=0.75,rely=0.9575,anchor=tk.CENTER)

button10 = tk.Button(frame2, text = 'Submit',
                      bg = 'DarkOrchid4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DarkOrchid4", 
                      highlightbackground="DarkOrchid4", 
                      borderwidth=4,fg='white',width=20,height=2,command=button10clicker)
button10.place(relx=0.25,rely=0.9575,anchor=tk.CENTER)

button11 = tk.Button(frame2, text = 'Export',
                      bg = 'DarkOrchid4',
                      bd =  15, 
                      highlightthickness=3, 
                      highlightcolor="DarkOrchid4", 
                      highlightbackground="DarkOrchid4", 
                      borderwidth=4,fg='white',width=20,height=1,command=button11clicker)
button11.place(relx=0.75,rely=0.5625,anchor=tk.CENTER)

#defining a function to indicate the number of results displayed
def label16updater():
    global tree,frame2
    count = 0
    for i in tree.get_children():
        count = count + 1
    displayresults = str(count) + " compound(s) displayed"
    label16 = tk.Label(master=frame2,text=displayresults,bg='indian red')
    label16.place(relx=0.25,rely=0.5625,anchor='center')
    frame2.after(1000, label16updater)

#calling the function to show number of compounds displayed
label16updater()

# Starting the Tkinter event loop
window.mainloop()

#closing the cursor and the connection to the database when the user closes the GUI.
cur1.close()
con1.close()