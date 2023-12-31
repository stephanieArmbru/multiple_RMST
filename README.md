\section*{R Code instruction}

The code pipeline computes the Restricted Mean Survival Time (RMST) for multiple events. 
In particular, it estimates the joint covariance while accounting for the correlation in time-to-event between multiple events for the combined RMST difference between treatment and control group. 

\subsection*{0. Set-up}
The code pipeline comes in a compressed folder \textit{Pipeline}, containing
\begin{itemize}
    \item \textit{Figures} folder: pipeline saves survival plots 
    \item \textit{Data} folder: pipeline saves truncation time, survival estimates and RMST calculations
    \item \textit{functions\textunderscore pipeline.R}: contains all required functions for the computation; no action required
    \item \textit{Pipeline.R}: main script which if run produces output; input required 
\end{itemize}

If run for the first time, R must install packages. 
\begin{center}
\begin{BVerbatim}
install.packages(c("tidyverse", "survRM2"))
devtools::install_github(repo = "zrmacc/MRMST", 
                         force = TRUE)
\end{BVerbatim}
\end{center}

Once the packages are installed, they need to be loaded every time the code is run. 
\begin{center}
\begin{BVerbatim}
library(tidyverse)
library(survRM2)
library(MRMST)
\end{BVerbatim}
\end{center}

In the script \textit{Pipeline.R}, the working directory should be set to the location in which the Pipeline folder is saved. 

\begin{center}
\begin{BVerbatim}
setwd("~/.../Pipeline")
\end{BVerbatim}
\end{center}



\subsection*{1. Data} 
It works with any data set as long as it contains, (1) \textbf{time-to-event} $T_i$ with the column name \verb|t2EVENTNAME| and (2) \textbf{status indicator} $C_i$ with the column name \verb|EVENTNAME| or \verb|c4EVENTNAME| which equals $C_i = 1$ if the event was observed and $C_i = 0$ if it was not.

The data cannot contain any missing values. 
If an event was not observed, its censoring time must be inserted and the count variable set to 0. 

The dataset can be saved in the \textit{Pipeline} folder. \\
The path to the dataset must be specified in 
\begin{center}
\begin{BVerbatim}
data_raw <- read_csv("path_to_file.csv", 
                     show_col_types = FALSE)
\end{BVerbatim}
\end{center}


\subsection*{2. Input}
The time-to-event as well as the status for the multiple events of interest must be inserted in 
\begin{center}
\begin{BVerbatim}
times_to_events <- c(...) \\
statuses_for_events <- c(...)
\end{BVerbatim}
\end{center}

The pipeline constructs composite endpoints from the events in question and some competing risk, e.g. death.
The time-to-event and status column name for the competing risk must be inserted in 
\begin{center}
\begin{BVerbatim}
time_to_death <- \\
status_for_death <- 
\end{BVerbatim}
\end{center}

The column name for the variable which indicates the arm or treatment group must be inserted in 
\begin{center}
\begin{BVerbatim}
arm <- 
\end{BVerbatim}
\end{center}

The variable \verb|alpha| can be adjusted to the desired confidence interval level. 
Its default is \verb|alpha = 0.05|. 

The KM survival curve contains the number of patients at risk for arm 1 and 0 for a subset of times across the observed trial duration. 
The density of this grid can be set by inserting in
\begin{center}
\begin{BVerbatim}
nar_grid <- 
\end{BVerbatim}
\end{center}
The default is \verb|nar_grid = 10|. 

The times at which the number at risk are reported are rounded. 
The rounding precision can be set by inserting in 
\begin{center}
\begin{BVerbatim}
round_to <-
\end{BVerbatim}
\end{center}
The default is \verb|round_to = 50|. 


\subsection*{3. Workflow}
After completing the input, the remaining R script can be run. 

It produces composite endpoints and estimates the Kaplan-Meier survival curves. 
The survival curves are drawn and saved in the \textit{Figures} folder.
They are also delivered as output by the code. 
To display the plots, the user must hit <RETURN> in the R console. 
The plots illustrate the KM survival estimate for each arm, the pointwise $1-\alpha$ confidence interval (i.e. for each time point individually constructed) and a table of patients at risk. 


The default truncation time is set to the maximum truncation time. 
The truncation time across all event types can be manually set, based on the survival curves, by substituting in line 106, 
\begin{center}
\begin{BVerbatim}
chosen_tau <- 
\end{BVerbatim}
\end{center}

The remaining code calculates the RMST difference across the multiple endpoints as well as construct confidence intervals and performs the WLW procedure \cite{wei1989regression,wei1997overview}.
The results are saved as RData files in the \textit{Data} folder. 

