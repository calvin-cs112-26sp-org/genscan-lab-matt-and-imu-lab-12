/* genScan.cpp searches a genome file F for a given subsequence SS,
 *  and outputs the number of occurrences of SS within F.
 *  This solution:
 *   reads F char-by-char, using ifstream::get() and string::operator+=(),
 *   and scans using string::substr() and string::operator==().
 *
 * @author: Joel Adams, Calvin University, June 2025
 *
 * Usage: ./genScan <genomeFile> <subsequence> 
 *
 * Precondition: genomeFile is a text file containing 
 *                an organism's DNA sequence in *plain* format.
 */

#include <iostream>            // cout, etc.
#include <fstream>             // ifstream, etc.
#include <cstdlib>             // exit()
#include <omp.h>               // omp_get_wtime()
//#include <algorithm>           // search()
//#include "OO_MPI_IO.h"         // ParallelReader
using namespace std;


/* retrieve inputs from command line
 * @param: argc, an int
 * @param: argv, a char**
 * @param: file, a string&
 * @param: subSeq, a string&
 * Precondition: argc and arg are the main function parameters
 *                 for retrieving command line values
 *            && the user has entered a file name and subsequence
 *                 on the command line
 * Postcondition: file contains the name of the file the user entered
 *            && subSeq contains the name of the subsequence the user entered
 *                (i.e., what they want to search for).
 */
void processCommandLineArgs(int argc, char** argv, 
                            string& file, string& subSeq) {
    if (argc == 3) {
        file = string(argv[1]);
        subSeq = string(argv[2]);
    } else {
        cerr << "\n *** Usage: ./genomeScanner <fileName> <subSequence>\n\n";
        exit(1);
   }
}

/* read a specified number of bytes from a file into a string.
 * @param: fileName, a const string&
 * @param: seq, a string&
 * Precondition: fileName contains the name of the input file
 *           && seq is the empty string into which the chars should be read.
 * Postcondition: seq contains the chars from the input file.
 */
void readFile(const string& fileName, string& seq) {
  ifstream fin(fileName.data()); // open file
  if (! fin.is_open()) {         // check it
     cerr << "\n *** Unable to open '" << fileName
          << "' as input file\n\n";
     exit(1);
  } 
  fin.seekg(0, std::ios::end);   // jump to end
  long N = fin.tellg();          // N = fin's offset 
  fin.seekg(0, std::ios::beg);   // undo jump
  seq.resize(N);                 // resize seq to N
  fin.get((char*)seq.data(), N); // read N bytes
  fin.close();                   // clean up
}

/* scan a string containing a genetic sequence for a subsequence
 * @param: seq, a string
 * @param: subSeq, a string
 * Precondition: seq contains a genetic sequence
 *           && subSeq contains a subsequence
 * Postcondition: the function returns the number of occurrences
             of subSeq within seq.
 */
long scan(const string& seq, const string& subSeq) {
   long skip = subSeq.size();
   long count = 0;
   size_t pos = seq.find(subSeq, 0);
   while (pos != string::npos) {
      ++count; 
      pos = seq.find(subSeq, pos+skip);
   } 
   return count;
}


/* output results of scan
 * @param: subSeq, a string
 * @param: numSubSeqs, a long
 * @param: inputTime, a double
 * @param: scanTime, a double
 * @param: totalTime, a double
 * Precondition: subSeq contains the user's subsequence
 *           && numSubSeqs == the number of subSeq occurrences in the sequence
 * Postcondition: subSeq, numSubSeqs, inputTime, scanTime, and totalTime
 *                have been displayed with appropriate labels.
 */   
void printResults(const string& subSeq, long numSubSeqs,
                   double inputTime, double scanTime, double totalTime) {
  cout << "\n1 thread found " << numSubSeqs << " occurrence"
        << ((numSubSeqs == 1) ? "" : "s")
        << " of '" << subSeq 
        << "'\n\tRead Time \tScan Time \tTotal Time\n"
        << fixed << '\t' << inputTime << '\t' << scanTime << '\t' 
        << totalTime << "\n\n";
}

// --------- main function ---------------------
int main(int argc, char** argv) { 
    string fileName;
    string subSeq;
    string dna;

    double startTotalTime = omp_get_wtime();
    processCommandLineArgs(argc, argv, fileName, subSeq);

    double startReadTime = omp_get_wtime();
    readFile(fileName, dna);
    double readTime = omp_get_wtime() - startReadTime;

    double startScanTime = omp_get_wtime();
    long count = scan(dna, subSeq);
    double scanTime = omp_get_wtime() - startScanTime;
    double totalTime = omp_get_wtime() - startTotalTime;

    //printResults(subSeq, count, 0.0, 0.0, 0.0);
    printResults(subSeq, count, readTime, scanTime, totalTime);
}

