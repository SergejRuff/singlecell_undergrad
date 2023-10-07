# 1 Uebung

#1 Öffnen Sie R Studio. Mit welcher R-Version arbeiten Sie?
R.version

#2 Erstellen Sie ein neues R-Script. Benennen Sie es mit "R_script_2022-10-05" und speichern Sie es als
#".R"-file.

#File-> save as....

#3 Beginnen Sie das Script mit einem Befehl, der automatisch den gesamten R-Workspace cleared. Tipp:
# Machen Sie sich mit den Befehlen rm() und ls() vertraut.

#Beispiel

data1<- c(1,2,3)
data2<-("Goodbye World. I will be deleted soon")

ls()     #ls shows list of objects in enviroment 

rm(list=ls()) #list= will specify what has to be removed.it treats the enviroment as a list and remove everything.
#way more efficient than typing every object in the ()

# rm() allows you to remove objects

rm(list=setdiff(ls(),"data2")) 

#explain: rm removes the objects in enviroment. setdiff()shows the difference between objects.
#Meaning everything except data2 will be selected for removel. See next line what setdiff does
setdiff(ls(),"data2") #the difference here when selecting for data2 is data1. data1 will be removed

#Alternativ

remove(data2) # oder rm(a)

#4 In welchem Working Directory / Arbeitspfad arbeiten Sie gerade?

getwd()

# um neuen Directory zu setzen

setwd() #New path in ()

###############################################################################

#Übung 2

# 1 Führen Sie folgende Rechnung in R aus: df
#  log((256 ??? 721 + 374/100)5 ??? 5)

log(((256 * 721) + ((374/100))*5) - 5)

#2 Implementieren Sie die Rechnung in R, indem Sie den folgenden Variablen die angegebenen Werte zuweisen.
#Anschließend weisen Sie das Ergebnis der Variable "ergebnis" zu.

a = 256 
b = 721 
c = 374 
d = 100
e = 5
ergebnis <- log(((256 * 721) + ((374/100))*5) - 5)

###############################################################################

#Übung 3

#1  Welcher Objekttyp ist Ihr Objekt "ergebnis"?

class(ergebnis)

#2 Was ergibt class(ergebnis)?

#numeric

# 3 Weisen Sie dem Objekt "name" Ihren Vornamen zu. Welcher Objekttyp ist "name"?

name<-"Sergej"
class(name)
#character/String

#Weisen Sie dem Objekt "eins" die Zahl 1 zu. Welcher Objekttyp ist "eins"? Ändern Sie den Objekttyp zu
#"character", indem Sie die Funktion "as.character()" verwenden. Überprüfen Sie den neuen Objekttyp und
#geben Sie "eins" aus. Was fällt Ihnen auf?

eins<-1  
class(eins)     # numeric
einsy<-as.character(1)
class(einsy)
einsy
# Antowrt: "", um zu zeigen,dass es sich hierbei um einen String handelt.

#4 Ändern Sie den Objekttyp von Objekt "eins" in "factor". Überprüfen Sie den neuen Objekttyp und geben
#Sie "eins" aus. Was fällt Ihnen auf? Geben Sie die Level von "eins" aus.

factor<-as.factor(eins)
class(factor)
factor
#levels:1

#Factors are variables in R which take on a limited number of different values, such as month names or 
#weekdays. That is useful in statistical analysis, as it allows to store and treat such data correctly and 
#use it in different models.
#For example, let's say our dataframe needs a gender column, which can take one of the following two 
#values: "Male" or "Female".
#We create a factor using the factor function, passing it a vector of values that we need in our column:
gender <- factor(c("Male", "Female")) 
gender
#Answer:
#[1] Male Female Male 
#Levels: Female Male


# 6 Überlegen Sie sich eine logische Verknüpfung (Tipp: ein Beispiel einer logischen Verknüpfung ist "1 < 10")
#und weisen Sie diese einem Objekt zu. Welche Klasse hat dieses Objekt und welchen Wert hat dieses Objekt?

Bool<-10<24
class(Bool)
Bool

# 7 Überlegen Sie, welches Ergebnis der nachfolgende Ausdruck bringt (TRUE/FALSE) und probieren Sie
#das Verhalten aus:

8.9999999999999999999999999 == 9.00000000000000000000000000000

#####################################################################################

# Übung 4

# 1 Erstellen Sie Vektor v. Welcher Objekttyp ist es

v<-c(1,2,6,10,4)
class(v)
#numeric
v[2]

# 2 4.2 Berechnen Sie die Summe über die ersten drei Vektorelemente von ???v. Berechnen Sie den Mittelwert und
#die Standardabweichung über alle Vektorelemente von ???v.

sum(v[1:3])
mean(v)
sd(v)
# alternativ
sqrt(sum((v-mean(v))^2/(length(v)-1)))

# 3 Multiplikation mit Skalar: Multiplizieren Sie ???v mit 10. Was ist das Ergebnis? Weisen Sie das Ergebnis
#einem neuen Objekt "w" zu.

w<- v*10
w

# 4 Multiplizieren Sie ???v mit ???w. Was ist das Ergebnis? Weisen Sie das Ergebnis einem neuen Objekt "z" zu.
#Beschreiben Sie den Rechenvorgang, wie er in R ausgeführt wurde.

z<-w*v
z

#5 erzeuge neuen Vector mit rep(). Multipliziere y und v

y<-c(rep(5,5))
y
y*v
