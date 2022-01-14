import re
import time
import datetime
import importlib
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from random import *
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from run_saver import Run_Saver
from _thread import start_new_thread

# matplotbib in tkinter (gui framework)
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from one_county_single_simulation_runner import single_simulation_run


from keras.models import load_model

# "main method"
# all gui elements and global used variables are in this function
def make_gui():
    # "highest" gui element ==> window
    global window
    window = Tk()
    window.resizable(False, False)
    window.title("Corona-Machine-Learning-Model Generator")
    window.geometry("400x720")



    ######################          Loading and Running         ################################


    button_path = Button(window, text="Select Dataset:", command=command_select_dataset)
    button_path.place(x=10, y=10, width=150, height=20)

    global text_path
    text_path = Text(window)
    text_path.place(x=170, y=10, width=200, height=23)
    text_path.insert("1.0", "data/epi_data_final.csv")


    label_eta = Label(window, anchor="w", text="Learningrate:")
    label_eta.place(x=10, y=40, width=150, height=20)

    global text_eta
    text_eta = Text(window)
    text_eta.place(x=170, y=40, width=200, height=20)
    text_eta.insert("1.0", "0.001")


    label_epochs = Label(window, anchor="w", text="Epochs:")
    label_epochs.place(x=10, y=70, width=150, height=20)

    global text_epochs
    text_epochs = Text(window)
    text_epochs.place(x=170, y=70, width=200, height=20)
    text_epochs.insert("1.0", "80")


    label_batch_size = Label(window, anchor="w", text="Batch size:")
    label_batch_size.place(x=10, y=100, width=150, height=20)

    global text_batch_size
    text_batch_size = Text(window)
    text_batch_size.place(x=170, y=100, width=200, height=20)
    text_batch_size.insert("1.0", "512")


    label_metrics = Label(window, anchor="w", text="Metrics:")
    label_metrics.place(x=10, y=130, width=150, height=20)

    global text_metrics
    text_metrics = Text(window)
    text_metrics.place(x=170, y=130, width=200, height=20)
    text_metrics.insert("1.0", "mae, mse")


    label_loss = Label(window, anchor="w", text="Lossfunction:")
    label_loss.place(x=10, y=160, width=150, height=20)

    global text_loss
    text_loss = Text(window)
    text_loss.place(x=170, y=160, width=200, height=20)
    text_loss.insert("1.0", "mae")


    label_name = Label(window, anchor="w", text="Run name:")
    label_name.place(x=10, y=190, width=150, height=20)

    global text_name
    text_name = Text(window)
    text_name.place(x=170, y=190, width=200, height=23)
    text_name.insert("1.0", "best_model")


    label_ml_model = Label(window, anchor="w", text="Choose Model:")
    label_ml_model.place(x=10, y=220, width=150, height=20)

    global combobox_choose_ml_model
    combobox_choose_ml_model = ttk.Combobox(window, values=[],  state="readonly")
    combobox_choose_ml_model.place(x=170,y=220, width=200, height=20)
    combobox_choose_ml_model.bind("<<ComboboxSelected>>", command_update_choose_ml_model_selected)


    global button_load
    button_load = Button(window, text="Load dataset", command=command_load)
    button_load.place(x=30,y=245, width=100, height=20)

    global button_run
    button_run = Button(window, text="Run", command=command_run, state=DISABLED, bg="red")
    button_run.place(x=150,y=245, width=100, height=20)



    ############################        Selection and view of selected data         ####################################


    label_runs = Label(window, anchor="w", text="Currently running: ")
    label_runs.place(x=10, y=280,height=25)

    global listbox_runs
    listbox_runs = Listbox(window)
    listbox_runs.bind("<<ListboxSelect>>", command_select_runs)
    scrollbar_runs = Scrollbar(listbox_runs)
    scrollbar_runs.pack(side=RIGHT, fill=Y)
    listbox_runs["yscrollcommand"] = scrollbar_runs.set
    listbox_runs.place(x=10,y=300, width=360, height=50)
    scrollbar_runs.config(command=listbox_runs.yview)




    label_imported_runs = Label(window, anchor="w", text="Finished runs: (Double-Click for evaluation)")
    label_imported_runs.place(x=10, y=370,height=25)

    global listbox_imported_runs
    listbox_imported_runs = Listbox(window)
    listbox_imported_runs.bind("<<ListboxSelect>>", command_select_imported_runs)
    listbox_imported_runs.bind("<Double-1>", command_select_and_evaluate_imported_runs)
    scrollbar_imported_runs = Scrollbar(listbox_imported_runs)
    scrollbar_imported_runs.pack(side=RIGHT, fill=Y)
    listbox_imported_runs["yscrollcommand"] = scrollbar_imported_runs.set
    listbox_imported_runs.place(x=10,y=390, width=360, height=80)
    scrollbar_imported_runs.config(command=listbox_imported_runs.yview)





    label_selected_name = Label(window, anchor="w", text= "Run name: ")
    label_selected_name.place(x=10,y=505,width=100,height=20)

    label_selected_dataset_name = Label(window, anchor="w", text= "Datasetname: ")
    label_selected_dataset_name.place(x=10,y=530,width=100,height=20)

    label_selected_n_epochs = Label(window, anchor="w", text= "Epochs: ")
    label_selected_n_epochs.place(x=10,y=560,width=100,height=20)

    label_selected_n_batch_size = Label(window, anchor="w", text= "Batch size: ")
    label_selected_n_batch_size.place(x=10,y=580,width=100,height=20)

    label_selected_eta = Label(window, anchor="w", text= "Eta: ")
    label_selected_eta.place(x=10,y=600,width=100,height=20)

    label_selected_loss_func = Label(window, anchor="w", text= "Lossfunction: ")
    label_selected_loss_func.place(x=10,y=620,width=100,height=20)

    label_selected_metrics_funcs = Label(window, anchor="w", text= "Metrics: ")
    label_selected_metrics_funcs.place(x=10,y=640,width=100,height=20)

    label_selected_running_status = Label(window, anchor="w", text= "Running status:")
    label_selected_running_status.place(x=10,y=660,width=110,height=20)

    label_selected_ml_model_name = Label(window, anchor="w", text = "ML-Model name:")
    label_selected_ml_model_name.place(x=10,y=680,width=110,height=20)


    global text_selected_name
    text_selected_name = Text(window)
    text_selected_name.place(x=170, y=505, width=200, height=23)

    global text_selected_dataset_name
    text_selected_dataset_name = Text(window)
    text_selected_dataset_name.place(x=170, y=530, width=200, height=23)

    global text_selected_n_epochs
    text_selected_n_epochs = Text(window)
    text_selected_n_epochs.place(x=170, y=560, width=200, height=20)

    global text_selected_n_batch_size
    text_selected_n_batch_size = Text(window)
    text_selected_n_batch_size.place(x=170, y=580, width=200, height=20)

    global text_selected_eta
    text_selected_eta = Text(window)
    text_selected_eta.place(x=170, y=600, width=200, height=20)

    global text_selected_loss_func
    text_selected_loss_func = Text(window)
    text_selected_loss_func.place(x=170, y=620, width=200, height=20)

    global text_selected_metrics_funcs
    text_selected_metrics_funcs = Text(window)
    text_selected_metrics_funcs.place(x=170, y=640, width=200, height=20)

    global text_selected_running_status
    text_selected_running_status = Text(window)
    text_selected_running_status.place(x=170, y=660, width=200, height=20)

    global text_selected_ml_model_name
    text_selected_ml_model_name = Text(window)
    text_selected_ml_model_name.place(x=170, y=680, width=200, height=20)


    #########################################       Evaluation          ########################################



    # frame_evaluation background
    fbg= "gray82"
    # frame where all evaluation elements are in
    global frame_evaluation
    frame_evaluation = Frame(window, bg=fbg)
    frame_evaluation.place(x=400,y=10,width=1090, height=700)


    label_evaluation_title = Label(frame_evaluation, text="Evaluation", font=("Arial",20), bg=fbg, anchor="center")
    label_evaluation_title.place(x=0,y=5,width=700, height=30)

    label_evaluation_calculation_time = Label(frame_evaluation, text="Calculation time: ", bg=fbg, anchor="w")
    label_evaluation_calculation_time.place(x=10, y=50, width=150, height=20)

    label_evaluation_choose_metrics = Label(frame_evaluation, text="Choose metrics: ", bg=fbg, anchor="w")
    label_evaluation_choose_metrics.place(x=10, y=70, width=130, height=20)

    global combobox_choose_metrics
    combobox_choose_metrics = ttk.Combobox(frame_evaluation, values=[],  state="readonly")
    combobox_choose_metrics.place(x=150,y=70, width=100, height=20)
    combobox_choose_metrics.bind("<<ComboboxSelected>>", command_update_choose_metrics_selected)


    global frame_metrics_plot
    frame_metrics_plot = Frame(frame_evaluation, bg=fbg)
    frame_metrics_plot.place(x=60, y=140, width=550,height=550)


    global label_evaluation_calculation_time_value
    label_evaluation_calculation_time_value = Label(frame_evaluation, text="", bg=fbg, anchor="w")
    label_evaluation_calculation_time_value.place(x=150, y=50, width=110,height=20)

    label_evaluation_show_metric_train = Label(frame_evaluation, bg=fbg, text="Train: ")
    label_evaluation_show_metric_train.place(x=10,y=95,width=100,height=20)

    label_evaluation_show_metric_test = Label(frame_evaluation, bg=fbg, text="Test: ")
    label_evaluation_show_metric_test.place(x=10,y=115,width=100,height=20)


    global label_evaluation_show_metric_train_value
    label_evaluation_show_metric_train_value = Label(frame_evaluation, bg=fbg, text="")
    label_evaluation_show_metric_train_value.place(x=100,y=95,width=140,height=20)

    global label_evaluation_show_metric_test_value
    label_evaluation_show_metric_test_value = Label(frame_evaluation, bg=fbg, text="")
    label_evaluation_show_metric_test_value.place(x=100,y=115,width=140,height=20)



    button_random_test = Button(frame_evaluation, text="Random Test", command=command_random_test)
    button_random_test.place(x=290, y=50, width=100,height=20)



    global combobox_choose_compartment
    combobox_choose_compartment = ttk.Combobox(frame_evaluation, values=[], state=DISABLED)
    combobox_choose_compartment.place(x=410,y=50,width=200,height=20)
    combobox_choose_compartment.bind("<<ComboboxSelected>>", command_update_choose_y_value)
    combobox_choose_compartment["values"] = ["Susceptible", "Exposed", "Carying", "Infected","Hospitalized", "ICU", "Recovered", "Deaths", "All", "All excepts susceptibles"]
    combobox_choose_compartment.current(9)

    global combobox_choose_age_group
    combobox_choose_age_group = ttk.Combobox(frame_evaluation, values=[], state=DISABLED)
    combobox_choose_age_group.place(x=410,y=75,width=200,height=20)
    combobox_choose_age_group.bind("<<ComboboxSelected>>", command_update_choose_y_value)
    combobox_choose_age_group["values"] = ["0-4", "5-14", "15-34", "35-59","60-79", "80+", "All"]
    combobox_choose_age_group.current(2)




    global label_comparisson
    label_comparisson = Label(frame_evaluation, bg=fbg, font=("Arial",15),text="Performance comparisson test")
    label_comparisson.place(x=640,y=50,width=450,height=40)

    global button_select_dampings
    button_select_dampings = Button(frame_evaluation, text="Select dampings:", command=command_select_dampings)
    button_select_dampings.place(x=640,y=140,width=150,height=20)

    global listbox_dampings
    listbox_dampings = Listbox(frame_evaluation)
    scrollbar_dampings = Scrollbar(listbox_dampings)
    scrollbar_dampings.pack(side=RIGHT, fill=Y)
    listbox_dampings["yscrollcommand"] = scrollbar_dampings.set
    listbox_dampings.place(x=790,y=140, width=300, height=50)
    scrollbar_dampings.config(command=listbox_dampings.yview)


    global label_choose_start_values
    label_choose_start_values = Label(frame_evaluation, bg=fbg, font=("Arial",12),text="Select example day from where the start values sould be used")
    label_choose_start_values.place(x=640,y=220,width=450,height=30)



    #here you can select the day where you want to start the simulation from
    global combobox_choose_start_values
    combobox_choose_start_values = ttk.Combobox(frame_evaluation, values=[i for i in range(180)])
    combobox_choose_start_values.place(x=640, y=250, width=450, height=30)
    combobox_choose_start_values.bind("<<ComboboxSelected>>", command_update_choose_start_values)
    combobox_choose_start_values.current(0)

    global button_compare
    button_compare = Button(frame_evaluation, text="Compare",command=command_compare)
    button_compare.place(x=640, y=300, width=100, height=20)


    #Labels to show results
    global label_show_time_numeric_result
    label_show_time_numeric_result = Label(frame_evaluation, bg=fbg, font=("Arial",12),text="Numeric simulation time: ", anchor="w")
    label_show_time_numeric_result.place(x=640,y=380,width=250,height=20)

    global label_time_numeric_result
    label_time_numeric_result = Label(frame_evaluation, bg=fbg, font=("Arial",12),text="-", anchor="w")
    label_time_numeric_result.place(x=900,y=380,width=200,height=20)



    global label_show_time_mlm_result
    label_show_time_mlm_result = Label(frame_evaluation, bg=fbg, font=("Arial",12),text="Machine learning prediction time: ", anchor="w")
    label_show_time_mlm_result.place(x=640,y=410,width=250,height=20)

    global label_time_mlm_result
    label_time_mlm_result = Label(frame_evaluation, bg=fbg, font=("Arial",12),text="-", anchor="w")
    label_time_mlm_result.place(x=900,y=410,width=200,height=20)


    ######################################      Global used variables       ###########################################


    global list_imported_runs
    list_imported_runs = []

    global list_runs
    list_runs = []

    # selected evaluating index
    global sidx
    sidx = -1

    global list_mlm
    list_mlm = [None]

    global list_mlm_imports
    list_mlm_imports = []
    #list_mlm_imports = [importlib.import_module("2DampingsLearning").__getattribute__("ML_Model")]

    #combobox_choose_ml_model["values"] = [str(list_mlm_imports[0]).replace(">","").replace("<","").replace("class","")]
    #combobox_choose_ml_model.current(0)
    global mlm_index
    mlm_index = 0


    ####################################        Menu options        ####################################################

    menu_main = Menu(window)

    menu_file = Menu(menu_main, tearoff=0)
    menu_file.add_command(label="Import Runs", command=command_import)
    menu_file.add_command(label="Add ML-Model", command=command_add_ml_model)

    global gpu_setter
    gpu_setter=BooleanVar()
    gpu_setter.set(False)
    global menu_options
    menu_options = Menu(menu_main, tearoff=0)
    menu_options.add_checkbutton(label="Deactivate GPU", onvalue=1, offvalue=0, variable=gpu_setter)

    menu_main.add_cascade(label="File", menu=menu_file)
    menu_main.add_cascade(label="Options", menu=menu_options)




    window.config(menu = menu_main)
    window.mainloop()


######################################################################################################################################################################



# import runs from directory
def command_import():
    import_file_path = filedialog.askopenfilenames(initialdir="./models/", title="Select run",filetypes=(("Run files", "*.run"),))
    for path in import_file_path:
        if(len(path) < 3):
            print("WARNING: Path to short")
            messagebox.showwarning("WARNING", "Import failed, one selected path is too short")
            return
        run = Run_Saver.load(path)
        if(run.getName() in [item.getName() for item in list_runs] or run.getName() in [item.getName() for item in list_imported_runs]):
            print("WARNING: Name already exists")
            messagebox.showwarning("WARNING", "Import failed, already imported a file with the same name")
        else:
            list_imported_runs.append(run)
            listbox_imported_runs.insert(END, run.getName())


def command_add_ml_model():
    ml_model_path = filedialog.askopenfilename(initialdir="./", title="Select Machine Learning Model",filetypes=(("Python files", "*.py"),))
    if(len(ml_model_path) < 4):
        return
    ml_model_path = re.sub(".*?/","",ml_model_path,flags=re.DOTALL).replace(".py","")

    list_mlm_imports.append(importlib.import_module(ml_model_path).__getattribute__("ML_Model"))
    list_mlm.append(None)
    combobox_choose_ml_model["values"] = [str(item).replace("<","").replace(">","").replace("class","") for item in list_mlm_imports]





# select dataset from directory
def command_select_dataset():
    dataset_path = filedialog.askopenfilename(initialdir="./data/", title="Select Dataset",filetypes=(("CSV files", "*.csv"),))
    text_path.delete("1.0","end")
    text_path.insert("1.0", dataset_path)



# load dataset in program
def command_load():
    button_load.configure(bg="red")
    file_name = text_path.get("1.0","end").replace("\n","")
    list_mlm[mlm_index] = list_mlm_imports[mlm_index](file_name)
    button_load.configure(bg="green")
    button_load["text"] = "Load again"
    button_run.configure(bg="green")
    button_run["state"] = NORMAL



# run ML-Model with given parameters
def command_run():
    button_run["state"] = DISABLED
    button_run.configure(bg="red")
    eta = float(text_eta.get("1.0","end").replace("\n",""))
    n_epochs = int(text_epochs.get("1.0","end").replace("\n",""))
    n_batch_size = int(text_batch_size.get("1.0","end").replace("\n",""))
    metrics_funcs = text_metrics.get("1.0","end").replace(" ","").replace("\n","").split(",")
    loss_func = text_loss.get("1.0","end").replace("\n","")
    dataset_name = text_path.get("1.0", "end").replace("\n","")
    ml_model_name = str(list_mlm_imports[mlm_index]).replace("<","").replace(">","").replace("class","")
    name = text_name.get("1.0","end").replace("\n","")
    name = "./models/" + name

    if(name in [item.getName() for item in list_runs] or name in [item.getName() for item in list_imported_runs]):
        print("WARNING: Name already exists")
        messagebox.showwarning("WARNING", "Name already exists, run doesnt start")
    else:
        menu_options.entryconfig("Deactivate GPU", state="disabled")
        start_new_thread(new_ML_Model_run, (list_mlm[mlm_index], n_epochs, n_batch_size, name, eta, loss_func, metrics_funcs, dataset_name, ml_model_name, gpu_setter))
        time.sleep(1)

    button_run.configure(bg="green")
    button_run["state"] = NORMAL



# select currently running run and show data
def command_select_runs(selection):
    run_name = listbox_runs.get(ANCHOR)
    for item in list_runs:
        if(item.getName() == run_name):
            if(item.getRunning_status() == 1):
                hideEvaluation()
            printInfos(item)
            break


def command_select_imported_runs(selection):
    # selection
    global sidx
    imported_run_name = listbox_imported_runs.get(ANCHOR)
    for item in list_imported_runs:
        if(item.getName() == imported_run_name):
            printInfos(item)
            sidx = list_imported_runs.index(item)
            break


# select finished run and show data
def command_select_and_evaluate_imported_runs(selection):
    # selection
    global sidx
    imported_run_name = listbox_imported_runs.get(ANCHOR)
    for item in list_imported_runs:
        if(item.getName() == imported_run_name):
            printInfos(item)
            sidx = list_imported_runs.index(item)
            break

    # evaluation 
    showEvaluation()
    label_evaluation_calculation_time_value["text"] = list_imported_runs[sidx].getCalculation_time()
    combobox_choose_metrics["values"] = list_imported_runs[sidx].getMetrics_funcs()
    combobox_choose_metrics.current(0)
    command_update_choose_metrics_selected(None)

    combobox_choose_compartment.current(9)
    combobox_choose_compartment["state"] = DISABLED

    combobox_choose_age_group.current(2)
    combobox_choose_age_group["state"] = DISABLED



# metrics value and plot update after new picked metrics
def command_update_choose_metrics_selected(selection):
    for widget in frame_metrics_plot.winfo_children():
        widget.destroy()

    label_evaluation_show_metric_train_value["text"] = round(list_imported_runs[sidx].getScores()[0][combobox_choose_metrics.current()],5)
    label_evaluation_show_metric_test_value["text"] = round(list_imported_runs[sidx].getScores()[1][combobox_choose_metrics.current()],5)

    figure_metrics = Figure(figsize=(3,3), dpi=100)
    a = figure_metrics.add_subplot(111)
    a.plot(list_imported_runs[sidx].getHistory()[combobox_choose_metrics.current()], color="blue", label="Train metric")
    a.set_xlabel("Epochs")
    a.set_ylabel("ERROR")
    a.legend()

    canvas = FigureCanvasTkAgg(figure_metrics, frame_metrics_plot)
    canvas.get_tk_widget().pack(side=TOP,fill=BOTH, expand=True)

    toolbar = NavigationToolbar2Tk(canvas, frame_metrics_plot)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP,fill=BOTH, expand=True)



def command_update_choose_ml_model_selected(selection):
    global mlm_index
    mlm_index = combobox_choose_ml_model.current()
    if(list_mlm[mlm_index] == None):
        button_load["text"] = "Load dataset"
        button_run["state"] = DISABLED
        button_run.configure(bg="red")
    else:
        button_load["text"] = "Load again"
        button_run["state"] = NORMAL
        button_run.configure(bg="green")


def command_update_choose_y_value(selection):
    showRandomTest()



def command_random_test():
    combobox_choose_compartment["state"] = "readonly"
    combobox_choose_age_group["state"] = "readonly"
    global random_test_idx
    random_test_idx = randint(0, len(list_imported_runs[sidx].getTestY())-1)
    showRandomTest()


def command_select_dampings():
    import_damping_paths = filedialog.askopenfilenames(initialdir="./predictions/dampings/", title="Select damping",filetypes=(("CSV files", "*.csv"),))
    listbox_dampings.select_clear
    if(len(import_damping_paths) == 0):
        messagebox.showwarning("WARNING", "Damping Import failed")
        return
    for damping in import_damping_paths:
        listbox_dampings.insert(END, damping.split("/")[len(damping.split("/"))-1])

def command_update_choose_start_values(selection):
    print("hi")


def command_compare():
    dampings = getDampings()
    #number of days which should be simulated. Depends on ML model
    num_sim_days = 10

    t1 = time.time()
    #make numeric simulation
    #simulation_result[0] = X, simulation_result[1] = Y
    simulation_result, sim_time = single_simulation_run(combobox_choose_start_values.current(), dampings, num_sim_days)
    numeric_sim_time = time.time() - t1


    t1 = time.time()
    #make prediciton with machine learning model
    #load existing model
    model = load_model(list_imported_runs[sidx].getName() + ".h5")
    #prepare cnn_input
    cnn_input = [np.array([[np.array(row) for row in dampings[i]]]) for i in range(0,len(dampings), 2)]
    #prepare mlp input
    mlp_input = [list([simulation_result[0][i] for i in range(1,len(dampings),2)]) + list(simulation_result[0][len(dampings):])] #for item in simulation_result[0]] #np.array([np.array([item[1] for item in dampings] + getStartValues())])
    #conversion to normalized mlp inputs
    mlp_input = list_imported_runs[sidx].getXMLPScaler().transform(mlp_input)
    #combine inputs
    x_input = [mlp_input] + cnn_input
    #make predictions
    prediction = model.predict(x_input)
    #transform normalited outputs to usefull outputs
    mlm_result = list_imported_runs[sidx].getYScaler().inverse_transform(prediction)


    mlm_time = time.time() - t1

    print(mlm_result)
    print("len")
    print(len(mlm_result[0]))
    print(len(simulation_result[1]))
    from sklearn.metrics import mean_absolute_error
    mae = mean_absolute_error(simulation_result[1], mlm_result[0])#1/len(mlm_result)* sum([((simulation_result[1][i] - mlm_result[0][i])/simulation_result[1][i])*100 for i in range(len(mlm_result[0]))])
    print("MAE: " + str(mae))

    label_time_mlm_result["text"] = mlm_time
    label_time_numeric_result["text"] = numeric_sim_time

########################################################################################################################################


def showRandomTest():
    global random_test_idx

    for widget in frame_metrics_plot.winfo_children():
        widget.destroy()

    combobox_compartment_idx = combobox_choose_compartment.current()
    combobox_age_group_idx = combobox_choose_age_group.current()

    # format: test_values[day][age_group][compartment]
    timesteps = int(len(list_imported_runs[sidx].getTestY()[random_test_idx])/(8*6))
    test_values = [[[ list_imported_runs[sidx].getTestY()[random_test_idx][compartment+(8*age_group)+(48*timestep)] for compartment in range(8)] for age_group in range(6)] for timestep in range(timesteps)]
    predicted_values = [[[ list_imported_runs[sidx].getPredtest()[random_test_idx][compartment+(8*age_group)+(48*timestep)] for compartment in range(8)] for age_group in range(6)] for timestep in range(timesteps)]


    figure_metrics = Figure(figsize=(3,3), dpi=100)
    a = figure_metrics.add_subplot(111)

    colors=["blue", "red", "green", "purple", "orange", "cornflowerblue", "orchid", "sienna", "black", "cyan", "chocolate", "indigo", "navy", "lightpink", "lightcoral", "brown"]





    if(combobox_age_group_idx < 6 and combobox_compartment_idx < 8):
        test_values_select = [item[combobox_age_group_idx][combobox_compartment_idx] for item in test_values]#[day[combobox_compartment_idx] for day in test_values[combobox_age_group_idx]]
        predicted_values_select = [item[combobox_age_group_idx][combobox_compartment_idx] for item in predicted_values]#[day[combobox_compartment_idx] for day in predicted_values[combobox_age_group_idx]]

        a.plot(test_values_select, color=colors[0], label="Real", dashes=[3,3])
        a.plot(predicted_values_select, color=colors[0], label="Predicted")

    elif(combobox_age_group_idx == 6 and combobox_compartment_idx < 8):
        test_values_select = [[day[i][combobox_compartment_idx] for day in test_values] for i in range(6)]
        predicted_values_select = [[day[i][combobox_compartment_idx] for day in predicted_values] for i in range(6)]

        for i in range(6):
            a.plot(test_values_select[i], color=colors[i], label=combobox_choose_age_group["values"][i]  + " R", dashes=[3,3])
            a.plot(predicted_values_select[i], color=colors[i], label=combobox_choose_age_group["values"][i]  + " P")

    elif(combobox_age_group_idx < 6 and combobox_compartment_idx == 8):
        test_values_select = [[day[combobox_age_group_idx][i] for day in test_values] for i in range(8)]
        predicted_values_select = [[day[combobox_age_group_idx][i] for day in predicted_values] for i in range(8)]

        for i in range(8):
            a.plot(test_values_select[i], color=colors[i], label=combobox_choose_compartment["values"][i]  + " R", dashes=[3,3])
            a.plot(predicted_values_select[i], color=colors[i], label=combobox_choose_compartment["values"][i]  + " P")

    elif(combobox_age_group_idx < 6 and combobox_compartment_idx == 9):
        test_values_select = [[day[combobox_age_group_idx][i] for day in test_values] for i in range(1,8,1)]
        predicted_values_select = [[day[combobox_age_group_idx][i] for day in predicted_values] for i in range(1,8,1)]

        for i in range(7):
            a.plot(test_values_select[i], color=colors[i], label=combobox_choose_compartment["values"][i+1]  + " R", dashes=[3,3])
            a.plot(predicted_values_select[i], color=colors[i], label=combobox_choose_compartment["values"][i+1]  + " P")
    elif(combobox_age_group_idx == 6 and combobox_compartment_idx >= 8):
        print("Not meaningful")
        return



    a.set_xlabel("Time")
    a.set_ylabel("Count of each compartment")
    a.legend()


    canvas = FigureCanvasTkAgg(figure_metrics, frame_metrics_plot)
    canvas.get_tk_widget().pack(side=TOP,fill=BOTH, expand=True)

    toolbar = NavigationToolbar2Tk(canvas, frame_metrics_plot)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP,fill=BOTH, expand=True)



# print run infos from selected run
def printInfos(run):
    text_selected_name.delete("1.0","end")
    text_selected_name.insert("1.0",run.getName())

    text_selected_dataset_name.delete("1.0","end")
    text_selected_dataset_name.insert("1.0",run.getDataset_name())

    text_selected_n_epochs.delete("1.0","end")
    text_selected_n_epochs.insert("1.0",run.getN_epochs())

    text_selected_n_batch_size.delete("1.0","end")
    text_selected_n_batch_size.insert("1.0",run.getN_batch_size())

    text_selected_eta.delete("1.0","end")
    text_selected_eta.insert("1.0",run.getEta())

    text_selected_loss_func.delete("1.0","end")
    text_selected_loss_func.insert("1.0",run.getLoss_func())

    text_selected_metrics_funcs.delete("1.0","end")
    text_selected_metrics_funcs.insert("1.0",run.getMetrics_funcs())

    text_selected_running_status.delete("1.0","end")
    text_selected_running_status.insert("1.0",run.getRunning_status())

    text_selected_ml_model_name.delete("1.0","end")
    text_selected_ml_model_name.insert("1.0",run.getMl_model_name())

# hide evaluation frame and reset all shown infos
def hideEvaluation():
    frame_evaluation.place_forget()
    combobox_choose_metrics["values"] = []
    combobox_choose_metrics.set("")
    for widget in frame_metrics_plot.winfo_children():
        widget.destroy()
    label_evaluation_show_metric_test_value["text"] = ""
    label_evaluation_show_metric_train_value["text"] = ""
    window.geometry("400x720")


# show evaluation frame 
def showEvaluation():
    window.geometry("1500x720")


# new run with ML-Model and given parameters
def new_ML_Model_run(mlm, n_epochs, n_batch_size, name, eta, loss_func, metrics_funcs, dataset_name, ml_model_name, gpu_setter):
    run = Run_Saver(n_epochs, n_batch_size, name, eta, loss_func, metrics_funcs, dataset_name, ml_model_name)
    list_runs.append(run)
    listbox_runs.insert(END,name)

    if(gpu_setter.get()):
        mlm.deactivateGPU()
    else:
        mlm.activateGPU()
    scores, history, predtest, calculation_time = mlm.run_Model(n_epochs=n_epochs, n_batch_size=n_batch_size, name=name, eta=eta, loss_func=loss_func, metrics_funcs=metrics_funcs)
    run.setResults(scores, history, predtest, mlm.getTestY(), mlm.getYScaler(), mlm.getXMLPScaler() ,calculation_time)

    run.save(name)
    list_runs.remove(run)

    index = -1
    for i in range(listbox_runs.size()):
        if(listbox_runs.get(i,i+1) == name):
            index = i

    listbox_runs.delete(index, index+1)

    list_imported_runs.append(run)
    listbox_imported_runs.insert(END, name)

def getDampings():
    dampings = []
    for item in list(listbox_dampings.get(0,END)):
        tmp = open("./predictions/dampings/"+item,"r").read()
        arr = [[float(tmp) for tmp in row.split(" ") if tmp!=""] for row in tmp.split("\n")]
        dampings.append(arr[1:7])
        dampings.append(arr[0][0])
    return dampings



make_gui()
