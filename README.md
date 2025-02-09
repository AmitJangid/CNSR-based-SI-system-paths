# CNSR-based-SI-system-paths

The role of promiscuous molecular recognition in the evolution of RNase-based self-incompatibility
Keren Erez (1), Amit Jangid (1), Ohad Noy Feldheim (2) and Tamar Friedlander (1)

	https://www.nature.com/articles/s41467-024-49163-7
	
    (1) The Robert H. Smith Institute of Plant Sciences and Genetics in Agriculture
        Faculty of Agriculture, The Hebrew University of Jerusalem,
        P.O. Box 12 Rehovot 7610001, Israel
    (2) The Einstein Institute of Mathematics, Faculty of Natural Sciences,
        The Hebrew University of Jerusalem, Jerusalem 9190401, Israel.
       
    Correspondence: tamar.friedlander@mail.huji.ac.il.

Keren Erez, Amit Jangid: Equal contrubution
###############################################################################################################

Source Code:
    The folder 'source_code' contains all the Python scripts which are required to regenerate the data, and figures in the manuscript.

Simulated dataset:
    The folder 'simulated_dataset' contains final output datasets and Python script to regenerate the figures in the main manuscript and also       
    in the supplementary.



The following are the requirements for the Python scripts:

1. System requirements
    a) All the scripts are written in Python and are simulated in Python 3.8
    b) Following are libraries followed by their versions used in the simulation in Python
        numpy 1.20.3
        pickle 4.0
        matplotlib 3.5.0
        seaborn 0.12.2
        jupyter 1.0.0
        os Miscellaneous operating system interfaces [https://docs.python.org/3/library/os.html#]
    c) The scripts were simulated on a local server with a clock speed of 2.4 GHz. It does not require any non-standard hardware.

2. Installation guide:
    a) Python and other libraries that are used to simulate the code are standard and can be installed in many different ways.
       One of the simplest ways is to install via ANACONDA [https://www.anaconda.com/].
    b) We did not use any specific software for simulation.

3. Demo:
    a) Here we provide the data and the Python scripts to regenerate the figures in the main manuscript in the folder 'simulated_dataset'. 
       Data and corresponding scripts are provided in the individual sub-folder named fig{i} respective to the figure number in the main manuscript.
    b) Running the scripts will provide the figures which are in the main manuscript and in the supplementary.
    c) The run time for only these scripts (for given datasets) is under 3 munites. 
       To regenerate the whole data for a single run, it takes 90 to 120 hrs on a clock speed of 2.4 GHz.

4. Instructions for use:
    a) No specific software has been used. We use custom codes which are explained above.
    b) We also provide a Jupyter notebook in the folder 'source code' which explains which Python script is needed to run and in which order.

###############################################################################################################  
