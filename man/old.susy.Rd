\name{old.susy}
\alias{old.susy}
\title{ Surrogate Synchrony }
\description{
  Cross-correlations are computed up to a specific lag in seconds \code{Maxlag}, then aggregated within a chosen \code{Segment} of e.g. 30s. Aggregation is performed by transforming correlations to Fisher's Z, then computing mean Z in each segment, then across all segments of the time series.
Segment shuffling is used to create surrogate time series, on which the same computations are run. This provides effect sizes \code{ES}.

\code{File}: The pairs of time series are in the columns of the file, variable names are in the header line. If \code{Automatic} is clicked, the synchrony is computed of all adjacent pairs of columns in the file. If \code{Automatic} is unclicked, you may choose the two columns to be analyzed for synchrony, and two plots are additionally prepared.
SUSY provides two different synchrony measures of each twin time series: mean Z and ES of mean Z; mean absolute_Z and ES of mean absolute_Z.
}
\usage{
old.susy(
  filename="MEA-topstage-C3-Beethoven.txt", separator=" ", missings=".", datahead=TRUE,
  spalte1=2, spalte2=5, epoche=30, fps=15, maxlag=3*fps, automatic=FALSE,
  plotid=3, plotid2=5, pseudosVereinfachen=FALSE, anzahlPseudosGesamt=500,
  zufallsdaten=FALSE, anzahlzufallsdaten=10000, korrOhneBetragAnzeigen=TRUE,
  permutation=FALSE
)
}
\arguments{
  \item{filename}{ Character scalar, name of a source file. Defaults to \code{"MEA-topstage-C3-Beethoven.txt"}. }
  \item{separator}{ Character scalar, source file data separator. Defaults to \code{" "}. }
  \item{missings}{ Character scalar, source file data missings marker. Defalts to \code{"."}. }
  \item{datahead}{ Logical scalar, tells if source data has a header row. Defaults to \code{TRUE}. }
  \item{spalte1}{ Numeric, ???hier bezeichnete Spalte 1 wird mit Spalte 2 verglichen (falls automatic = FALSE)??. Defaults to \code{2}. }
  \item{spalte2}{ Numeric, ???. Defaults to \code{5}. }
  \item{epoche}{ Epochengroesse (SEGMENT SIZE) in s (Epoche darf nicht kleiner sein als Zeitfenster/fps)
Epoche*fps darf nicht groesser sein als die Anzahl Daten (timepoints) geteilt durch 2
Beispiel 1: timepoints: 1000, fps = 10, also darf Epoche nicht groesser sein als 50s (=1000/(10*2))
Beispiel 2: maxlag : 90, fps = 10, also darf Epoche nicht KLEINER sein als 9s (=90/10). }
  \item{fps}{ frames per second (SAMPLING RATE). }
  \item{maxlag}{ Zeitfenster (WINDOW SIZE) fuer CCF - > maximaler lag in Sekunden (nur ganze Zahlen). }
  \item{automatic}{ automatische Iteration? TRUE, FALSE
es wird die erste mit der zweiten Spalte verglichen, dann die dritte mit der vierten usw...
ein Spaltenpaar muss dieselbe Anzahl Daten haben, aber verschiedene Paare nicht.
Beachte: Es werden keine Plots ausgegeben falls automatic = TRUE. }
  \item{plotid}{ welche Art plot soll angezeigt werden (man kann den plot auch manuell ausfuehren (unten))
0 nichts, 1 GMcrosscorrs, 2 synchrony by epoch, 3 GM-Z, 4 Zeitreihenplot, 5 z not abs
wenn plotid2 == 0, dann kein zweiter plot }
  \item{plotid2}{ ... }
  \item{pseudosVereinfachen}{ Berechnung der pseudo-Daten vereinfachen (TRUE, FALSE)
Falls FALSE, anzahlPseudosGesamt irrelevant
anzahlPseudos: Anzahl der Surrogats-tests aller Epochen }
  \item{anzahlPseudosGesamt}{ ... }
  \item{zufallsdaten}{ Zufallsdaten (Normalverteilte Zahlen im Bereich [0,50] }
  \item{anzahlzufallsdaten}{ ... }
  \item{korrOhneBetragAnzeigen}{ Korrelationen ohne Betrag auch anzeigen? (FALSE/TRUE) }
  \item{permutation}{ Permuation: Synchronie f?r n(n-1)/2 Dyaden berechnen }
}
\value{
  \code{NULL} is returned. Plots are created as a side effect.
}
\examples{
if (file.exists("MEA-topstage-C3-Beethoven.txt"))
  old.susy()
}
\keyword{ data }
