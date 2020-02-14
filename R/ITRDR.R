getCount <- function(x)
{
   require(RCurl)
   require(XML)
   require(httr)
   x <- RCurl::curlEscape(x)
   #Sys.sleep(0.01)
   baseURL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi?retmode=xml&term="
   x <- paste(baseURL,x,sep="")
   n <- 0
   while(n!=10)
   {
       y <- try(GET(x,timeout(30)),silent=TRUE)
       y <- try(content(y,encoding="UTF-8"),silent=TRUE)
       y <- try(xmlParse(y),silent=TRUE)
       if (class(y)[1L]=="try-error") 
       {
           cat("ERROR1: ", y, "\n")
           Sys.sleep(10)
           print("reconnecting...")
           n <- n+1
           print(n)
       } 
       else 
       {
          break
       } 
   }
   y <- xmlToList(y)
   y <- y$eGQueryResult$ResultItem$Count
   as.numeric(y)
}