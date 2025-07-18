```

██╗  ██╗ ██████╗ ██╗      █████╗      █████╗ ███╗   ███╗██╗ ██████╗  ██████╗ ███████╗
██║  ██║██╔═══██╗██║     ██╔══██╗    ██╔══██╗████╗ ████║██║██╔════╝ ██╔═══██╗██╔════╝
███████║██║   ██║██║     ███████║    ███████║██╔████╔██║██║██║  ███╗██║   ██║███████╗
██╔══██║██║   ██║██║     ██╔══██║    ██╔══██║██║╚██╔╝██║██║██║   ██║██║   ██║╚════██║
██║  ██║╚██████╔╝███████╗██║  ██║    ██║  ██║██║ ╚═╝ ██║██║╚██████╔╝╚██████╔╝███████║
╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝    ╚═╝  ╚═╝╚═╝     ╚═╝╚═╝ ╚═════╝  ╚═════╝ ╚══════╝
                                                                                     

```


Hi there, remember me? Yeah, the guy who interned from May–July '25?  
Yes, it's *my* directory. Yes, my intern is over, and now you need to tune parameters yourself.  
Any doubts? No worries — this is my mail ID: **arnavabhishek166@gmail.com**, just ping me.

Meanwhile, you can check out my GitHub: [github.com/ArnavAbhishek](https://github.com/ArnavAbhishek),  
and do connect with me on LinkedIn: [linkedin.com/arnavabhishek](https://linkedin.com/arnavabhishek).  

It was nice working with you. :)

**Signing off...**  
Abhishek Arnav  
IIT Delhi

---

## 📂 Directory Structure

This is the structure of the directory you're currently in:

Final_working_directory/
├── oneparameter/
│   ├── pso.py
│   ├── assets/
│   │   ├── error.py
│   │   ├── model.py
│   │   ├── data.py
│   │   ├── file_completion.py
│   │   ├── abaqus_completion.py
│   │   └── objective.py
│   └── umat_directory/
│       ├── umat.for
│       └── ... (other UMAT-related files)
│
├── twoparameter/
│   ├── pso.py
│   ├── assets/
│   │   ├── error.py
│   │   ├── model.py
│   │   ├── data.py
│   │   ├── file_completion.py
│   │   ├── abaqus_completion.py
│   │   └── objective.py
│   └── umat_directory/
│       ├── umat.for
│       └── ... (other UMAT-related files)


---

###  To run any of the tuning simulations, just follow these steps:

> **DISCLAIMER**: The present code tries to tune **XTACUC1** and **CRSS** values.  
> If you wish to change the parameters, follow Step III.

---

####  STEP I

- Open CMD prompt in the `Final_working_directory`
- Paste this command: `"cd oneparameter"` or `"cd twoparameter"` depending on what you want to tune
- Open `pso.py` in the same directory, or type: `code pso.py` on the CMD prompt
- Verify the final input/umat files on:
  - lines 25 and 27 for one parameter
  - lines 27 and 31 for two parameters  
  (if you wish to change them, follow Step II)
- Change the real parameters to the true experimental values (lines 25–27 or 27–31)
- There is a segment called `initial_guess` — you may try changing that too
- You may also change hyperparameters of the PSO algorithm — I have tuned them, but feel free to improve them 
- Now paste this "abaqus python pso.py" command on the CMD prompt (You may add `interactive` at the end — I prefer data logging rather than filling buffer cache)

 The tuning algorithm will start and you can monitor progress in real time via `output.txt`.

---

####  STEP II

Want to tune a different `.inp` file?

- Remove the current `.inp` and paste the new one
- Open `pso.py` and modify **line 24** (rename the file)
- If you have a different `umat.for`, paste it in `umat_directory` and similarly update **line 28**

---

####  STEP III

Since you're here, you want to change something. If yes, don’t...  
Just kidding.

- Create a Git session (or at least a backup)
- Create a new branch and then modify files in the `assets` folder
- You’ll mostly need to change `error.py` or `data.py`

If you're tuning **other parameters**, follow the steps below:

---

### A) Modifying `data.py`

- Please be very careful — this took me **2 weeks** to write
- Look for the `extractData()` function — it's defined for **peak load extraction**
- You can define your own function to extract the quantity you wish to tune
- For `twoparameter`, I have defined **4 functions**, each with specific roles — feel free to repurpose them

---

### B) Modifying `error.py`

- Also took me **2 weeks**, because I was editing the wrong file for days 
- Look for the `modify()` function
- Change the files where you want to tune parameters — for instance:
- `"kMaterialParam.f"` at **line 152** for **XTACUC1**

---

All the best tuning your simulations!  
Let me know if you break something 
