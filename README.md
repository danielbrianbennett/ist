# ist
Code to calculate P(correct) in the Information Sampling Task, using MATLAB

These functions use MATLAB's built-in xlsread and xlswrite functions to import data from a raw CANTAB output file (the 'input_filename' argument) and export it in processed form (the 'output_filename' argument). A couple of notes: 
  (1) Both input and output files must be Microsoft Excel spreadsheets. 
  (2) Because of limitations in the MATLAB functions that handle reading and writing from excel files, this function will only work in Windows, and not on Macs. We wrote and tested it using MATLAB version R2015b and Microsoft Excel 2010 for Windows. While we would not expect performance to be substantially different in other versions of these programs, we have not tested these and make no guarantees regarding the software.

To help users test this code, we have provided an example of a raw CANTAB output file ('raw_test_data.xlsx'), which we have converted into a processed datafile called 'processed_test_data.xlsx' using the FormatCANTAB.m function. Users should replicate this on their own machines to check that the code is working as it should before applying it to other raw files.

Any questions please contact Dan Bennett: danielbrianbennett (at) gmail (dot) com

This code is distributed under a GNU GPL license as documented at http://www.gnu.org/licenses/gpl.txt. Software is distributed as is, completely without warranty or service support. The authors of the software are not liable for the condition or performance of the code, and hereby disclaim all implied warranties.
