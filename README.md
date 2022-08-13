# GNAR_models_for_Irish_COVID-19_data
R code for Dissertation "Generalised Network Time Series models: comparison of network construction approaches in modelling the COVID-19 incidence in Ireland" by Stephanie Armbruster


The structure of the code is modular. 
The folder \textbf{data} has two subfolders \textbf{COVID} and \textbf{Mumps} with all relevant datasets. \\
Folder \textbf{COVID} includes: 
\begin{itemize}
    \item COVID-19\_HPSC\_full\_week.csv: data set with daily cumulative COVID-19 incidence from \cite{COVID_data}
    \item ireland\_covid\_weekly.csv: data set after pre-processing
    \item ireland\_shapefile.shp: shapefile for plotting network maps, generated by saving data downloaded from the GADM database\footnote{available at \url{https://gadm.org}, accessed: 2022-07-27} 
    \item county\_towns\_ireland.csv: data set with coordinates for Irish counties from the Ireland Cities Database on simplemaps\footnote{available at \url{https://simplemaps.com/data/ie-cities}, accessed: 2022-07-10}
\end{itemize}
Folder \textbf{Mumps} includes: 
\begin{itemize}
    \item mumps.RData: data set with daily mumps incidence for the UK in 2005
    \item england\_wales\_shapefile.shp: shapefile for plotting network maps, generated by saving data downloaded from the GADM database\footnote{available at \url{https://gadm.org}, accessed: 2022-07-27} 
\end{itemize}

The file \textit{functions.R} includes all functions written to pre-process data, fit \code{GNAR} models, analyse residuals, plot BIC values etc. 
In order to incorporate the alternative weighting schemes, the functions \code{GNARdesign()}, \code{GNARfit()} and \code{NofNeighbours()} from the package \code{GNAR} had to be expanded. 
The adapted functions \code{GNARdesign\_weighting()}, \code{GNARfit\_weighting()} and \code{NofNeighbours\_named()} are included in the files \textit{GNARdesign.R}, \textit{GNARfit.R} and \textit{NofNeighbours.R} in the folder \textbf{GNAR}.

The file \textit{Replication\_Knight\_2016.R} replicates the results from \cite{knight2016modelling}. 

The file \textit{COVID\_Ireland\_data\_processing.R} performs some initial data exploration, the data aggregation to a weekly level and the necessary smoothing for extreme peaks in COVID-19 incidence.

The COVID-19 networks are constructed and the \code{GNAR} models are fit and analysed in the file \textit{COVID\_Ireland.R}. 

The file \textit{COVID\_Ireland\_regulations.R} explores \code{GNAR} models fitted to the 5 data subsets which are constructed by splitting the original data set according to COVID-19 regulations. 

