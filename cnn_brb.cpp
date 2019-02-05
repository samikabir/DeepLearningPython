#include <iostream>
#include <math.h> 
#include <fstream>
#include <string>

using namespace std;

float a1 = 0;
float a2 = 0;

float PM_H = 500.4;
float PM_M = 35.5;
float PM_L = 0.0;

float AQI_H = 500.0;
float AQI_M = 101.0;
float AQI_L = 0.0;

float H1 = 0.0;
float M1 = 0.0;
float L1 = 0.0; 

float H2 = 0.0;
float M2 = 0.0;
float L2 = 0.0;

int numberOfAntAttributes = 2;
float matchingDegree[9];
float relativeWeight = 1.0;
float totalWeight = 0;
float consequentBeliefDegree[27];
float updatedConsequentBeliefDegree[27];
float beliefDegreeChangeLevel = 0; 
float activationWeight[9];
float ruleWiseBeliefDegreeSum[9]; 
string line;
string cnn_mild;
string cnn_nominal;
string cnn_severe;
int counter = 0;
float normalized_cnn_severe_degree = 1;
float normalized_cnn_mild_degree = 1;
float normalized_cnn_nominal_degree = 1;
float cnn_pm25 = 1;
float aggregatedBeliefDegreeH = 1;
float aggregatedBeliefDegreeM = 1;
float aggregatedBeliefDegreeL = 1;
float finalAggregatedBeliefDegreeH = 1.0;
float finalAggregatedBeliefDegreeM = 1.0;
float finalAggregatedBeliefDegreeL = 1.0;
float brbH = 0;
float brbM = 0;
float brbL = 0;
int aqi = 1;

void transformInput1(float i);
void transformInput2(float i);

void ruleBase()
{
    consequentBeliefDegree[0] = 1;
    consequentBeliefDegree[1] = 0;
    consequentBeliefDegree[2] = 0;
    consequentBeliefDegree[3] = 0;
    consequentBeliefDegree[4] = 0.40;
    consequentBeliefDegree[5] = 0.60;
    consequentBeliefDegree[6] = 0;
    consequentBeliefDegree[7] = 0;
    consequentBeliefDegree[8] = 1;
    consequentBeliefDegree[9] = 1;
    consequentBeliefDegree[10] = 0;
    consequentBeliefDegree[11] = 0;
    consequentBeliefDegree[12] = 0;
    consequentBeliefDegree[13] = 1;
    consequentBeliefDegree[14] = 0;
    consequentBeliefDegree[15] = 0;
    consequentBeliefDegree[16] = 0.40;
    consequentBeliefDegree[17] = 0.60;
    consequentBeliefDegree[18] = 0;
    consequentBeliefDegree[19] = 0;
    consequentBeliefDegree[20] = 1;
    consequentBeliefDegree[21] = 0.60;
    consequentBeliefDegree[22] = 0.40;
    consequentBeliefDegree[23] = 0;
    consequentBeliefDegree[24] = 0;
    consequentBeliefDegree[25] = 0;
    consequentBeliefDegree[26] = 1;
}

void takeInput()
{
    
    cout<<"Insert value for PM2.5 (between 0 and 500.4 µg/m3): ";
    cin>>a1;
    
    //cout<<"Insert value for Air Quality Index - AQI (between 0 and 500): ";
    //cin>>a2;
    
    transformInput1(a1);
    //transformInput2(a2);
}

void transformInput1(float i)
{
    if (i >= PM_H)
    {
        H1 = 1;
        M1 = 0;
        L1 = 0;
    }
    else if (i == PM_M)
    {
        H1 = 0;
        M1 = 1;
        L1 = 0;
    }
    else if (i <= PM_L)
    {
        H1 = 0;
        M1 = 0;
        L1 = 1;
    }    
    else if ((i <= PM_H) && (i >= PM_M))
    {
        M1 = (PM_H-i)/(PM_H-PM_M);
        H1 = 1 - M1;
        L1 = 0.0; 
    }
    else if ((i <= PM_M) && (i >= PM_L))
    {
        L1 = (PM_M-i)/(PM_M-PM_L);
        M1 = 1 - L1; 
        H1 = 0.0; 
    }
}

void transformInput2(float i)
{
    if (i >= AQI_H)
    {
        H2 = 1;
        M2 = 0;
        L2 = 0;
    }
    else if (i == AQI_M)
    {
        H2 = 0;
        M2 = 1;
        L2 = 0;
    }
    else if (i <= AQI_L)
    {
        H2 = 0;
        M2 = 0;
        L2 = 1;
    }        
    else if ((i <= AQI_H) && (i >= AQI_M))
    {
        M2 = (AQI_H-i)/(AQI_H-AQI_M);
        H2 = 1 - M2;
        L2 = 0.0; 
    }
    else if ((i <= AQI_M) && (i >= AQI_L))
    {
        L2 = (AQI_M-i)/(AQI_M-AQI_L);
        M2 = 1 - L2; 
        H2 = 0.0; 
    }
}


void showTransformedInput()
{
    cout<< endl << "Transformed Input is as follow." << endl;
    cout<< "PM2.5 = {(H, " << H1 << "); (M, " << M1 << "); (L, " << L1 << ")}" << endl;
    //cout<< "AQI = {(H, " << H2 << "); (M, " << M2 << "); (L, " << L2 << ")}" << endl;
}

void calculateMatchingDegree()
{
    int increment = 0;
    float ti1[3] = {H1, M1, L1};
    float ti2[3] = {H2, M2, L2};
    
    for (int c = 0; c < 3; c++)
        for (int d = 0; d < 3; d++){ 
                //weight[increment] = ti1[c] * ti2[d] * ti3[e];
                matchingDegree[increment] = pow(ti1[c], relativeWeight) * pow(ti2[d], relativeWeight);
                increment++; 
            }
}
  
void calculateMatchingDegreeBrbCnn()
{ 
    int increment = 0;
    float ti1[3] = {H1, M1, L1};    
    float ti2[3] = {normalized_cnn_severe_degree, normalized_cnn_mild_degree, normalized_cnn_nominal_degree};
    
    for (int c = 0; c < 3; c++)
        for (int d = 0; d < 3; d++){ 
                //weight[increment] = ti1[c] * ti2[d] * ti3[e];
                matchingDegree[increment] = pow(ti1[c], relativeWeight) * pow(ti2[d], relativeWeight);
                increment++; 
            }
}

void showMatchingDegree()
{
    int track = 1; 
    //cout << endl << "Matching degrees of the rules are as follow." << endl; 
    for (int counter = 0; counter < 9; counter++)
    {
        //cout<< "Matching Degree of Rule " << track << " = " << matchingDegree[counter] << endl;
        track++; 
    }
}

void showActivationWeight()
{
    int trace = 1; 
    for (int x = 0; x < 9; x++)
    {
        totalWeight += matchingDegree[x];
    }
    
    //cout << endl << "Activation Weights of the rules are as follow."<< endl; 
    
    for (int counter = 0; counter < 9; counter++)
    {
        activationWeight[counter] = matchingDegree[counter]/totalWeight; 
        //cout<< "Activation weight of Rule " << trace << " = " << activationWeight[counter] << endl;
        trace++;  
    }
}

void updateBeliefDegree()
{
    int update = 0;
    float sumAntAttr1 = 1;
    float sumAntAttr2 = 1;    
    
    if ((H1 + M1 + L1) < 1)
    {
        sumAntAttr1 = H1 + M1 + L1;
        update = 1;
    }
    
    if ((H2 + M2 + L2) < 1)
    {
        sumAntAttr2 = H2 + M2 + L2;
        update = 1;
    }
    
    if (update == 1)
    {
        beliefDegreeChangeLevel = (sumAntAttr1 + sumAntAttr2)/numberOfAntAttributes;
        //cout << "Belief Degree Change Level = " << beliefDegreeChangeLevel << endl;
        for (int go = 0; go < 27; go++)
        {
            consequentBeliefDegree[go] = beliefDegreeChangeLevel * consequentBeliefDegree[go];
            //cout << "Updated Consequent Belief Degree : " << consequentBeliefDegree[go] << endl;
        }
    }
    else
    {
        //cout << endl << "No upgradation of belief degree required." << endl;
    }
    
}

void aggregateER()
{
    int parse = 0;
    int move1 = 0; 
    int move2 = 1; 
    int move3 = 2; 
    int action1 = 0;
    int action2 = 1;
    int action3 = 2;
    
    float part11 = 1;
    float part12 = 1;
    float part13 = 1;
    float part1 = 1;
    float part2 = 1;
    float value = 1;
    float meu = 1;
    
    float numeratorH1 = 1;
    float numeratorH2 = 1;
    float numeratorH = 1;
    float denominatorH1 = 1;
    float denominatorH = 1;
    
    float numeratorM1 = 1;
    float numeratorM = 1;
    
    float numeratorL1 = 1;
    float numeratorL = 1;
    
    float utilityScoreH = 1;
    float utilityScoreM = 0.5;
    float utilityScoreL = 0;
    float crispValue = 1;
    float degreeOfIncompleteness = 1;
    float utilityMax = 1;
    float utilityMin = 1;
    float utilityAvg = 1; 
    
    for (int t = 0; t < 9; t++)
    {
        parse = t * 3;
        ruleWiseBeliefDegreeSum[t] = consequentBeliefDegree[parse] + consequentBeliefDegree[parse+1] + consequentBeliefDegree[parse+2];
    }
    
    for (int rule = 0; rule < 9; rule++){
        part11 *= (activationWeight[rule] * consequentBeliefDegree[move1] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        move1 += 3;
    }
     
    for (int rule = 0; rule < 9; rule++){
        part12 *= (activationWeight[rule] * consequentBeliefDegree[move2] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        move2 += 3;
    }

    for (int rule = 0; rule < 9; rule++){
        part13 *= (activationWeight[rule] * consequentBeliefDegree[move3] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        move3 += 3;
    }
    
    part1 = (part11 + part12 + part13);
    
    for (int rule = 0; rule < 9; rule++){
        part2 *= (1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
    }    
    
    value = part1 - part2;
    
    meu = 1/value;

    for (int rule = 0; rule < 9; rule++){
        numeratorH1 *= (activationWeight[rule] * consequentBeliefDegree[action1] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        action1 += 3;
    }
    
    for (int rule = 0; rule < 9; rule++){
        numeratorH2 *= (1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
    }    
    
    numeratorH = meu * (numeratorH1 - numeratorH2);
    
    for (int rule = 0; rule < 9; rule++){
        denominatorH1 *= (1 - activationWeight[rule]);        
    }
    
    denominatorH = 1 - (meu * denominatorH1);
    
    aggregatedBeliefDegreeH = (numeratorH/denominatorH);
    //cout << endl << "ER Aggregated Belief Degree for Severe Pollution: " << aggregatedBeliefDegreeH << endl;
    
    for (int rule = 0; rule < 9; rule++){
        numeratorM1 *= (activationWeight[rule] * consequentBeliefDegree[action2] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        action2 += 3;
    }
    
    numeratorM = meu * (numeratorM1 - numeratorH2); 
    aggregatedBeliefDegreeM = (numeratorM/denominatorH); 
    //cout << "ER Aggregated Belief Degree for Mild Pollution: " << aggregatedBeliefDegreeM << endl;
    
    for (int rule = 0; rule < 9; rule++){
        numeratorL1 *= (activationWeight[rule] * consequentBeliefDegree[action3] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        action3 += 3;
    }
     
    numeratorL = meu * (numeratorL1 - numeratorH2);
    aggregatedBeliefDegreeL = (numeratorL/denominatorH); 
    //cout << "ER Aggregated Belief Degree for Nominal Pollution: " << aggregatedBeliefDegreeL << endl;    
    
    if ((aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL) == 1){
        crispValue = (aggregatedBeliefDegreeH * utilityScoreH) + (aggregatedBeliefDegreeM * utilityScoreM) + (aggregatedBeliefDegreeL * utilityScoreL);
        //cout << "Crisp or numerical value is: " << crispValue << endl;        
        brbH = aggregatedBeliefDegreeH;
        brbM = aggregatedBeliefDegreeM;
        brbL = aggregatedBeliefDegreeL;    
        
        cout << endl << "ER Aggregated Belief Degree for Severe Pollution: " << aggregatedBeliefDegreeH << endl;
        cout << "ER Aggregated Belief Degree for Mild Pollution: " << aggregatedBeliefDegreeM << endl; 
        cout << "ER Aggregated Belief Degree for Nominal Pollution: " << aggregatedBeliefDegreeL << endl;
        //cout << "brbH: " << brbH << " brbM: " << brbM << " brbL: " << brbL <<endl;
    }

        
    else{
        
        degreeOfIncompleteness = 1 - (aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
        //cout << "Usassigned Degree of Belief: " << degreeOfIncompleteness << endl; 
        
        utilityMax = ((aggregatedBeliefDegreeH + degreeOfIncompleteness) * utilityScoreH + (aggregatedBeliefDegreeM*utilityScoreM) + (aggregatedBeliefDegreeL*utilityScoreL));
        
        utilityMin = (aggregatedBeliefDegreeH*utilityScoreH) + (aggregatedBeliefDegreeM*utilityScoreM) + (aggregatedBeliefDegreeL + degreeOfIncompleteness) * utilityScoreL;
        
        utilityAvg = (utilityMax + utilityMin)/2;
        
        //cout << "Maximum expected utility: " << utilityMax << endl;
        //cout << "Minimum expected utility: " << utilityMin << endl;  
        //cout << "Average expected utility: " << utilityAvg << endl; 
        
        cout << endl << "BRB ER Aggregated Belief Degrees considering degree of Incompleteness are as follow." << endl;
        
        finalAggregatedBeliefDegreeH = aggregatedBeliefDegreeH/(aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
        
        finalAggregatedBeliefDegreeM = aggregatedBeliefDegreeM/(aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
        
        finalAggregatedBeliefDegreeL = aggregatedBeliefDegreeL/(aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
            
        cout << endl << "BRB-generated Belief Degree for Severe Pollution: " << finalAggregatedBeliefDegreeH << endl; 
        cout << "BRB-generated Belief Degree for Mild Pollution: " << finalAggregatedBeliefDegreeM << endl; 
        cout << "BRB-generated Belief Degree for Nominal Pollution: " << finalAggregatedBeliefDegreeL << endl;  
        
        brbH = finalAggregatedBeliefDegreeH;
        brbM = finalAggregatedBeliefDegreeM;
        brbL = finalAggregatedBeliefDegreeL;  
        
        //cout << "brbH: " << brbH << " brbM: " << brbM << " brbL: " << brbL <<endl;
    } 
    
}

void takeCnnOutput() 
{ 
        
  ifstream myfile ("cnn_prediction.txt"); //cnn output 
  //ifstream myfile ("cnn_prediction1.txt"); //severe 408      
  //ifstream myfile ("cnn_prediction2.txt"); //nominal 36   
  //ifstream myfile ("cnn_prediction3.txt"); //mild 117        
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      //cout << line << '\n';
      if(counter == 0)
      {
         cnn_mild = line;        
      }
      else if (counter == 1)
      {
          cnn_nominal = line;        
      }
      else
          cnn_severe = line;        
        
      counter++;     
     
    }
    myfile.close();
  }
 
  else cout << "Unable to open file";  
  float a = stof(cnn_mild);
  float b = stof(cnn_nominal);
  float c = stof(cnn_severe);    
    
  float mild_degree = a/100;
  float nominal_degree = b/100;
  float severe_degree = c/100;
    
  float sum_degree = severe_degree + mild_degree + nominal_degree;
  
  normalized_cnn_severe_degree = severe_degree/sum_degree;
  normalized_cnn_mild_degree = mild_degree/sum_degree;   
  normalized_cnn_nominal_degree = nominal_degree/sum_degree;       
  
    
  if ((normalized_cnn_severe_degree > normalized_cnn_mild_degree) && (normalized_cnn_severe_degree > normalized_cnn_nominal_degree))
  {
      cnn_pm25 = (150.5 + 349.9*normalized_cnn_severe_degree) + ((150.4*normalized_cnn_mild_degree)/2);
      cout << endl << "PM2.5 computed by CNN: " << cnn_pm25 << " µg/m3 " << endl;          
  }
  else if ((normalized_cnn_nominal_degree > normalized_cnn_mild_degree) && (normalized_cnn_nominal_degree > normalized_cnn_severe_degree))
  { 
      cnn_pm25 = (35.4*(1 - normalized_cnn_nominal_degree)) + ((150.4*normalized_cnn_mild_degree)/2);            
      cout << endl << "PM2.5 computed by CNN: " << cnn_pm25 << " µg/m3 " << endl;          
  }
  else if ((normalized_cnn_mild_degree > normalized_cnn_severe_degree) && (normalized_cnn_mild_degree > normalized_cnn_nominal_degree))
  {
      if (normalized_cnn_severe_degree > normalized_cnn_nominal_degree)
      { 
        cnn_pm25 = (35.5 + 114.9*normalized_cnn_mild_degree) + ((500.4*normalized_cnn_severe_degree)/2);
        cout << endl << "PM2.5 computed by CNN: " << cnn_pm25 << " µg/m3 " << endl;          
      }
      
      else if ((normalized_cnn_nominal_degree > normalized_cnn_severe_degree))
      {
        cnn_pm25 = (35.5 + 114.9*normalized_cnn_mild_degree) + ((35.4*normalized_cnn_nominal_degree)/2);    
        cout << endl << "PM2.5 computed by CNN: " << cnn_pm25 << " µg/m3 " << endl;           
      }
  }
    
  cout << endl << "CNN-generated Belief Degree for Severe Pollution: " << normalized_cnn_severe_degree << endl;       
  cout << "CNN-generated Belief Degree for Mild Pollution: " << normalized_cnn_mild_degree << endl;   
  cout << "CNN-generated Belief Degree for Nominal Pollution: " << normalized_cnn_nominal_degree << endl;    
}

void aggregateER_BrbCnn()
{
    int parse = 0;
    int move1 = 0; 
    int move2 = 1; 
    int move3 = 2; 
    int action1 = 0;
    int action2 = 1;
    int action3 = 2;
    
    float part11 = 1;
    float part12 = 1;
    float part13 = 1;
    float part1 = 1;
    float part2 = 1;
    float value = 1;
    float meu = 1;
    
    float numeratorH1 = 1;
    float numeratorH2 = 1;
    float numeratorH = 1;
    float denominatorH1 = 1;
    float denominatorH = 1;
    
    float numeratorM1 = 1;
    float numeratorM = 1;
    
    float numeratorL1 = 1;
    float numeratorL = 1;
    
    float utilityScoreH = 1;
    float utilityScoreM = 0.5;
    float utilityScoreL = 0;
    float crispValue = 1;
    float degreeOfIncompleteness = 1;
    float utilityMax = 1;
    float utilityMin = 1;
    float utilityAvg = 1; 
    
    for (int t = 0; t < 9; t++)
    {
        parse = t * 3;
        ruleWiseBeliefDegreeSum[t] = consequentBeliefDegree[parse] + consequentBeliefDegree[parse+1] + consequentBeliefDegree[parse+2];
    }
    
    for (int rule = 0; rule < 9; rule++){
        part11 *= (activationWeight[rule] * consequentBeliefDegree[move1] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        move1 += 3;
    }
     
    for (int rule = 0; rule < 9; rule++){
        part12 *= (activationWeight[rule] * consequentBeliefDegree[move2] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        move2 += 3;
    }

    for (int rule = 0; rule < 9; rule++){
        part13 *= (activationWeight[rule] * consequentBeliefDegree[move3] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        move3 += 3;
    }
    
    part1 = (part11 + part12 + part13);
    
    for (int rule = 0; rule < 9; rule++){
        part2 *= (1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
    }    
    
    value = part1 - part2;
    
    meu = 1/value;

    for (int rule = 0; rule < 9; rule++){
        numeratorH1 *= (activationWeight[rule] * consequentBeliefDegree[action1] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        action1 += 3;
    }
    
    for (int rule = 0; rule < 9; rule++){
        numeratorH2 *= (1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
    }    
    
    numeratorH = meu * (numeratorH1 - numeratorH2);
    
    for (int rule = 0; rule < 9; rule++){
        denominatorH1 *= (1 - activationWeight[rule]);        
    }
    
    denominatorH = 1 - (meu * denominatorH1);
    
    aggregatedBeliefDegreeH = (numeratorH/denominatorH);
    //cout << endl << "ER Aggregated Belief Degree for Severe Pollution: " << aggregatedBeliefDegreeH << endl;
    
    for (int rule = 0; rule < 9; rule++){
        numeratorM1 *= (activationWeight[rule] * consequentBeliefDegree[action2] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        action2 += 3;
    }
    
    numeratorM = meu * (numeratorM1 - numeratorH2); 
    aggregatedBeliefDegreeM = (numeratorM/denominatorH); 
    //cout << "ER Aggregated Belief Degree for Mild Pollution: " << aggregatedBeliefDegreeM << endl;
    
    for (int rule = 0; rule < 9; rule++){
        numeratorL1 *= (activationWeight[rule] * consequentBeliefDegree[action3] + 1 - (activationWeight[rule] * ruleWiseBeliefDegreeSum[rule]));        
        action3 += 3;
    }
     
    numeratorL = meu * (numeratorL1 - numeratorH2);
    aggregatedBeliefDegreeL = (numeratorL/denominatorH); 
    //cout << "ER Aggregated Belief Degree for Nominal Pollution: " << aggregatedBeliefDegreeL << endl;    
    
    if ((aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL) == 1){
        crispValue = (aggregatedBeliefDegreeH * utilityScoreH) + (aggregatedBeliefDegreeM * utilityScoreM) + (aggregatedBeliefDegreeL * utilityScoreL);
        //cout << "Crisp or numerical value is: " << crispValue << endl;        
        brbH = aggregatedBeliefDegreeH;
        brbM = aggregatedBeliefDegreeM;
        brbL = aggregatedBeliefDegreeL;    
        
        cout << endl << "BRB-CNN integrated Belief Degree for Hazardous AQI: " << aggregatedBeliefDegreeH << endl;
        cout << "BRB-CNN integrated Belief Degree for Unhealthy AQI: " << aggregatedBeliefDegreeM << endl; 
        cout << "BRB-CNN integrated Belief Degree for Good AQI: " << aggregatedBeliefDegreeL << endl;
        //cout << "brbH: " << brbH << " brbM: " << brbM << " brbL: " << brbL <<endl; 
    }   
 
        
    else{
        
        degreeOfIncompleteness = 1 - (aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
        //cout << "Usassigned Degree of Belief: " << degreeOfIncompleteness << endl; 
        
        utilityMax = ((aggregatedBeliefDegreeH + degreeOfIncompleteness) * utilityScoreH + (aggregatedBeliefDegreeM*utilityScoreM) + (aggregatedBeliefDegreeL*utilityScoreL));
        
        utilityMin = (aggregatedBeliefDegreeH*utilityScoreH) + (aggregatedBeliefDegreeM*utilityScoreM) + (aggregatedBeliefDegreeL + degreeOfIncompleteness) * utilityScoreL;
        
        utilityAvg = (utilityMax + utilityMin)/2;
        
        //cout << "Maximum expected utility: " << utilityMax << endl;
        //cout << "Minimum expected utility: " << utilityMin << endl; 
        //cout << "Average expected utility: " << utilityAvg << endl; 
        
        cout << endl << "BRB-CNN integrated Belief Degrees considering degree of Incompleteness:" << endl;
        
        finalAggregatedBeliefDegreeH = aggregatedBeliefDegreeH/(aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
         
        finalAggregatedBeliefDegreeM = aggregatedBeliefDegreeM/(aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
        
        finalAggregatedBeliefDegreeL = aggregatedBeliefDegreeL/(aggregatedBeliefDegreeH + aggregatedBeliefDegreeM + aggregatedBeliefDegreeL);
                
        
        cout << endl << "BRB-CNN integrated Belief Degree for Hazardous AQI: " << finalAggregatedBeliefDegreeH << endl; 
        cout << "BRB-CNN integrated Belief Degree for Unhealthy AQI: " << finalAggregatedBeliefDegreeM << endl; 
        cout << "BRB-CNN integrated Belief Degree for Good AQI: " << finalAggregatedBeliefDegreeL << endl << endl; 
         
        brbH = finalAggregatedBeliefDegreeH;
        brbM = finalAggregatedBeliefDegreeM;
        brbL = finalAggregatedBeliefDegreeL;  
         
        //cout << "brbH: " << brbH << " brbM: " << brbM << " brbL: " << brbL <<endl;
        if ((finalAggregatedBeliefDegreeH > finalAggregatedBeliefDegreeM) && (finalAggregatedBeliefDegreeH > finalAggregatedBeliefDegreeL))
        {
            aqi = (201 + 299*finalAggregatedBeliefDegreeH) + ((200*finalAggregatedBeliefDegreeM)/2);
            cout << endl << "AQI predicted by BRB-CNN: " << aqi << endl;          
        }
        else if ((finalAggregatedBeliefDegreeL > finalAggregatedBeliefDegreeM) && (finalAggregatedBeliefDegreeL > finalAggregatedBeliefDegreeH))
        { 
            aqi = (100*(1 - finalAggregatedBeliefDegreeL)) + ((200*finalAggregatedBeliefDegreeM)/2);            
            cout << endl << "AQI predicted by BRB-CNN: " << aqi << endl;          
        }
        else if ((finalAggregatedBeliefDegreeM > finalAggregatedBeliefDegreeH) && (finalAggregatedBeliefDegreeM > finalAggregatedBeliefDegreeL))
        {
            if (finalAggregatedBeliefDegreeH > finalAggregatedBeliefDegreeL)
            { 
                aqi = (101 + 99*finalAggregatedBeliefDegreeM) + ((500*finalAggregatedBeliefDegreeH)/2);
                cout << endl << "AQI predicted by BRB-CNN: " << aqi << endl;           
            }
      
            else if ((finalAggregatedBeliefDegreeL > finalAggregatedBeliefDegreeH))
            { 
                aqi = (101 + 99*finalAggregatedBeliefDegreeM) + ((100*finalAggregatedBeliefDegreeL)/2);    
                cout << endl << "AQI predicted by BRB-CNN: " << aqi << endl;           
            }
        }
    }     
}

 
int main()
{
    ruleBase();
    takeInput();  
    showTransformedInput(); 
    /*calculateMatchingDegree();
    showMatchingDegree();
    showActivationWeight();
    updateBeliefDegree();
    aggregateER();*/
    
    takeCnnOutput(); 
    calculateMatchingDegreeBrbCnn();
    showMatchingDegree();
    showActivationWeight();
    updateBeliefDegree(); 
    aggregateER_BrbCnn(); //modify this function to calculate AQI. Do it on 28.11.2018 (Wed) morning 
     
    return 0;
}
