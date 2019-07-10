/* Author: Sucheendra kumar  NUS */ 

#include<stdio.h>
#include<string.h>
#include<time.h>
#include<stdlib.h>
#include <conio.h>
#define MAX 35
#define TIMEPOINT 100

int k=0;
int p=0,Formula_count=0,prop_count=0;
char label[35][250],label2[35][250],label3[250],proposition[5][15],num_char[4];
void parser(char in_str[][250]);
char *str_replace(char *str, char *orig, char *rep);
int cstring_cmp(const void *a, const void *b);
int substring(char *,char *);
char** unique_copy(char in_str[][250], int number);
char** a;
float brown[201][100][5];
int num_prop=0;

//array of structures for each time point.
typedef struct {
	char Formula[250];
	int Truth_flag;
	int evaluated;
	}LTL_Node_i[MAX];


//another struct to represent the nodes
typedef struct{
	LTL_Node_i LTL_Node_Element;
	int time_point;
	struct LTL_Node *nextptr;	
	struct LTL_Node *previousptr;
    }LTL_Node;

LTL_Node *getnode(int);/*For creating New Nodes(Malloc)*/
void new_node();
void update_node_Prop();
void update_node_remaining();
LTL_Node *head=NULL; // Initially the linked list is empty!
LTL_Node *tail=NULL;


int main()
{
    
	
	
	char inputexp[500]="(Px52=0,p>0.8)&F((Px52=4,p>0.8)&F(G(Px52=0,p>0.8)))";
	
	//oscillation cycle:	
	//Axin protein 1)"F((Px08=0,p>0.6)&F((Px08=2,p>0.6)&F((Px08=0,p>0.6)&F((Px08=2,p>0.6)&F((Px08=0,p>0.6)&F((Px08=2,p>0.6)&F((Px08=0,p>0.6)&F((Px08=2,p>0.6)&F((Px08=0,p>0.6)&F((Px08=2,p>0.6)&F(Px08=0,p>0.6)))))))))))";
	

	//brown model:
	//"(Px21=0,p>0.7)&F((Px21=3,p>0.7)&F(G(Px21=2,p>0.6)))";
	//"G(Px01=4,p>0.9)&G(Px02=4,p>0.9)";
	//"(Px27=0,p>0.7)&F(G(Px27=4,p>0.7))";

	//large Pathway model
	//"(Px22=0,p>0.7)&F((Px22=4,p>0.7)&F(G(Px22=0,p>0.7)))";
	//"(Px52=0,p>0.7)&F((Px52=4,p>0.7)&F(G(Px52=0,p>0.7)))";
	//"(Px57=0,p>0.7)&F((Px57=4,p>0.7)&F(G(Px57=0,p>0.7)))";
	//"(Px105=2,p>0.7)&F((Px105=4,p>0.7)&F(G(Px105=3,p>0.6)))";
	//"((Px52=4,p<0.1)U((Px22=4,p>0.8)&F(Px52=4,p>0.7)))"


	FILE *file1;
    int filet,filex,filev,j=0;
	int gone_in=0;
    float filep;
    int g,h;
	clock_t start,final;

	printf("**************************************************\n");
	printf("***** ------------------------------------ *******\n");
	printf("***  DBN Based Probabilistic Model Checker  ******\n");
	printf("***** ------------------------------------ *******\n");
	printf("**************************************************\n");
	printf("\n");
	printf("\n");
	printf("The Atomic Proposition should be of the form \n\n(Px[var]=[Interval],p[</>][probability value])\n");
	printf("\n");

	/* The user gets the LTL Expression to be checked */
	//printf("enter the LTL Formula that you want to evaluate:\n");
	//printf("\n");
	//scanf("%s",inputexp);
	//printf("\n");
	//printf("\n");



	start = clock();
	//file1= fopen("C:\\Documents and Settings\\Administrator\\My Documents\\Dropbox\\DBNmodelchecker\\m201_M_gpu.csv","r");
	//file1= fopen("C:\\Documents and Settings\\Administrator\\My Documents\\Dropbox\\DBNmodelchecker\\brown_M_gpu.csv","r");
	file1= fopen("C:\\Documents and Settings\\Administrator\\My Documents\\Dropbox\\DBNmodelchecker\\m88_M_gpu-0830.csv","r");
	//file1= fopen("C:\\Documents and Settings\\Administrator\\My Documents\\Dropbox\\DBNmodelchecker\\m88_M_gpu-3M.csv","r");
	//file1= fopen("C:\\Documents and Settings\\Administrator\\My Documents\\Dropbox\\DBNmodelchecker\\m88_M_cpu-30k.csv","r");

	/* read the contents of the file into a 3D Array */
	while(fscanf(file1,"%d,%d,%d,%f",&filet,&filex,&filev,&filep)== 4)
		{
                   brown[filex][filet][filev]=filep;
		}
     
     /* Close the file which contains the values of probabilities stored */
	fclose(file1);

	printf("The Model has %d Variables\n",filex);
	printf("Each Variable has %d intervals (0-%d)\n",filev+1, filev);
	printf("\n");
	printf("\n");

	for (g=0;g<35;g++)
	{
		memset(label[g],'\0',250);

	}


	printf("-----------------------------------------------\n");


	for(g=0;g<(int)strlen(inputexp);g++)
	{
		if(inputexp[g]=='G' || inputexp[g]=='F' || inputexp[g]=='X' || inputexp[g]=='&' || inputexp[g]=='|' || inputexp[g]=='P' || inputexp[g]=='!' || inputexp[g]=='U' || inputexp[g]=='-' )
		{
            Formula_count++;
		}
	}
	
    // Copy the entered first element into the array of substrings
    strcpy(label[0],inputexp);
            
	/* Start the clock to keep an account of the time taken */	
	//start = clock();
        
    // call the function parser that evaluates the subexpressions and stores them in an array of strings.
    parser(label);
	 
      
	for(g=0;g<Formula_count;g++)
	{
			strcpy(label2[g],label[g]);
	}

    // getting only the propositions in an array!
    for(g=0;g<Formula_count;g++)
	{
		if(label[g][0]=='P' && strlen(label[g])<15)
		{
			gone_in=0;
			for(h=0;h<15;h++)
			{
				if(strcmp(label[g],proposition[h])==0)
				{
					gone_in=1;
				}
			}
			if(gone_in==0)
			{
			strcpy(proposition[j],label[g]);
			num_prop++;
			j++;
			}
		}
	}
    
	//replacing the label array such with formulas containing numbers instead of propositions.
	for(g=0;g<Formula_count;g++)
	{
		j=0;
		while(strlen(proposition[j])!= 0)
		{
			if(strstr(label[g],proposition[j])!=NULL)
			{   
				itoa(j+1,num_char, 10);
				strcpy(label2[g],str_replace(label2[g],proposition[j],num_char));						
			}
			j++;
		}
	}

   qsort((void*)label2, Formula_count,250,cstring_cmp);

   a=unique_copy(label2,Formula_count);
    
   for(prop_count=-1;prop_count<TIMEPOINT-1;)
   {
    new_node();
   }

   update_node_Prop();
   update_node_remaining();

   if(head->LTL_Node_Element[0].Truth_flag==1)
   {
	   printf("The Formula %s is evaluated to be true\n",inputexp);
	   printf("-----------------------------------------------\n");
   }
   else if (head->LTL_Node_Element[0].Truth_flag==0)
   {
	   printf("The Formula %s is evaluated to be False\n",inputexp);
	   printf("-----------------------------------------------\n");
   }
   final=clock();
   printf("Time Taken : %f seconds",(double)(final-start)/CLOCKS_PER_SEC );
   getch();
   return 0;
          
  }  // End of Main     



 
/*******************************************************************
* NAME :            parser(char in_str[][25])
*
* DESCRIPTION :     Calculates all the substrings of the LTL string
*                   and stores them in an array of  strings.
*
* INPUTS :
*       PARAMETERS:
*           char     in_str[][25]          The array of substrings 
*      
* OUTPUTS :
*       RETURN : void
*           
* PROCESS :  1) accepts the LTL formula from the  user.
*            2) calculates all the substrings recursively.
*            3) The Array of strings has all the substrings.
*                  
*
* CHANGES :
* REF NO    DATE         WHO     DETAIL
*           30thjune2009 suchee  Original Code
            09thjuly2009 suchee  modified Code without recursion!
* *********************************************************************/

void parser(char in_str[][250])
    {
           
		int len,openbrac_count=0,len_counter[250]={0},no_brac=0;
		int i,from,some_count3,some_count4=0,and_count=0,and_count2=0,flag_count;
		char in_temp[250];

        while(in_str[k][0]!= '\0')
			{
				no_brac=0;
				memset(in_temp,'\0',250);
				
				while( no_brac!=1)
				{
				   len=strlen(in_str[k]) ;
				   for (i = 0; i < len; i++) 
				   {
                    len_counter[i] = 0;
                   }
				   openbrac_count=0;
				   i=0  ;
				   flag_count=0;
				   and_count2=0;
				   and_count=0;
				   some_count3=0;
				   some_count4=0;
				   from=0 ;
				   while(i<len)
							 {
								if(in_str[k][i]=='(')
								{
									for(openbrac_count=i;openbrac_count<len;openbrac_count++)
									{
										len_counter[openbrac_count]++;
									}
								}
								else if(in_str[k][i]==')')
								{
									for(openbrac_count=i;openbrac_count<len;openbrac_count++)
									{
										len_counter[openbrac_count]--;
									}
								}
								i++;
							}
				        if(in_str[k][0]=='(' && in_str[k][len-1]==')' && len_counter[0]==len_counter[len-2])
							{
									for(openbrac_count=0;openbrac_count<len-1;openbrac_count++)
											{
												if(len_counter[openbrac_count]<len_counter[0])
												{
													no_brac=1;
												}

											}
							}
						else if(in_str[k][0]!='(')
						{
                            no_brac=1;
						}
										
						   if(no_brac!=1)
						   {
							strncpy(in_temp,&in_str[k][1],strlen(in_str[k])-2); 
							memset(in_str[k],'\0',250);
							strcpy(in_str[k],in_temp);
						   }
		}




				    for(some_count3=0;some_count3<(int)strlen(in_str[k]);some_count3++)
					{
						if(in_str[k][some_count3]== '&' || in_str[k][some_count3]== '|' || in_str[k][some_count3]== 'U'|| in_str[k][some_count3]=='-')
							{
								and_count=and_count+1;
						    }
					}

				   for(some_count3=0;some_count3<len;some_count3++)
				   {
					   if(in_str[k][some_count3]=='&' || in_str[k][some_count3]=='|' || in_str[k][some_count3]== 'U')
					   {
						   and_count2=and_count2+1;
						   if(len_counter[some_count3]==0)
						   {
							   p=p+1;
							   strncpy(in_str[p],&in_str[k][some_count4],some_count3-some_count4); 
							   some_count4= some_count3+1;
							   flag_count=1;
						  	}
						   if(and_count2==and_count && flag_count==1)
								{
								  p=p+1;
								  strncpy(in_str[p],&in_str[k][some_count4],strlen(in_str[k])-some_count4);
								}
						 }
					   if(in_str[k][some_count3]=='-')
					   {
                          and_count2=and_count2+1;
						   if(in_str[k][some_count3-1]!='<' && in_str[k][some_count3+1]=='>')
						  {
						   if(len_counter[some_count3]==0)
							   {
								   p=p+1;
								   strncpy(in_str[p],&in_str[k][some_count4],some_count3-some_count4); 
								   some_count4= some_count3+2;
								   flag_count=1;
						  		}
						   if(and_count2==and_count && flag_count==1)
								{
									 p=p+1;
									 strncpy(in_str[p],&in_str[k][some_count4],strlen(in_str[k])-some_count4);
								}
						  }//end of in_str[k][some_count3-1]!='<' && in_str[k][some_count3-1]!='>'

						  else if(in_str[k][some_count3-1]=='<' && in_str[k][some_count3+1]=='>')
						  {
						   if(len_counter[some_count3]==0)
							   {
								   p=p+1;
								   strncpy(in_str[p],&in_str[k][some_count4],some_count3-some_count4-1); 
								   some_count4= some_count3+2;
								   flag_count=1;
						  		}
						   if(and_count2==and_count && flag_count==1)
								{
									 p=p+1;
									 strncpy(in_str[p],&in_str[k][some_count4],strlen(in_str[k]));
								}
						  }//end of in_str[k][some_count3-1]=='<' && in_str[k][some_count3-1]!='>'

						 }

					
				   }
			     if(flag_count!=1)
				 {
                 if (in_str[k][0]=='G'||in_str[k][0]=='X'||in_str[k][0]=='F')
					{
						p=p+1;
					    strncpy(in_str[p],&in_str[k][2],strlen(in_str[k])-3); 
					 }
				 if (in_str[k][0]=='!')
					{
					if(in_str[k][1]=='(')
							{
								p=p+1;
								strncpy(in_str[p],&in_str[k][2],strlen(in_str[k])-3); 
															   
							}
					else
							{
								p=p+1;
								strncpy(in_str[p],&in_str[k][1],strlen(in_str[k])); 
							}
														   
					}
				 
				 }
				 k++;								 
               }// End of while(in_str[k][0]!= '\0')    
           }//End of function parser

// replace function used to replace the propositions with numbers.
// char *str_replace(char *str, char *orig, char *rep)
// {
//   static char buffer[4096];
//   char *p;
// 
//   if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
//     return str;
// 
//   while(strstr(str, orig)!=NULL)
//   {
//   strncpy(buffer, str, p-str); // Copy characters from 'str' start to 'orig' st$
//   buffer[p-str] = '\0';
//   sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));
//   if(strstr(buffer, orig)!=NULL)
//   {
//   str_replace(buffer,orig,rep);
//   }
//   strcpy(str,buffer);
//   }
// 
//   return buffer;
//   }

		   char *str_replace(char *str, char *orig, char *rep)
		   {
			   static char buffer[4096];
			   char *p;
			   if(!(p = strstr(str, orig)))
				   return str;
			   while(strstr(str, orig)!=NULL)
			   {
				   p = strstr(str, orig);
				   strncpy(buffer, str, p-str);
				   buffer[p-str] = '\0';
				   sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));
				   strcpy(str,buffer);
			   }
			   return buffer;
		   }

//function to allocate nodes
LTL_Node *getnode(int prop_count)
{
	int g;
	LTL_Node *node;
	node=(LTL_Node *)malloc(sizeof(LTL_Node));
	node->time_point=prop_count;
	for(g=0;g<Formula_count;g++)
	{
	strcpy(node->LTL_Node_Element[g].Formula,a[g]);
	node->LTL_Node_Element[g].Truth_flag=0;
    node->LTL_Node_Element[g].evaluated=0;
	}
	for(g=Formula_count;g<MAX;g++)
	{
		memset(node->LTL_Node_Element[g].Formula,'\0',250);
	}
	node->previousptr=NULL;
    node->nextptr=NULL;

	return node;
}

void new_node()
{
	LTL_Node *node_1;

	if(head==NULL&&tail==NULL)
	{ 
		prop_count++;
		node_1=getnode(prop_count);
		head=node_1;
		tail=node_1;
		node_1->previousptr=NULL;
		node_1->nextptr=NULL;
	}
	else
	{ 
		prop_count++;
		node_1=getnode(prop_count);
		(LTL_Node*)node_1->previousptr=tail;
		(LTL_Node*)tail->nextptr=node_1;
		node_1->nextptr=NULL;
		tail=node_1;
	}
}

//updates all the propositions in the list as evaluated or not
void update_node_Prop()
{
   LTL_Node *node_2;
   int t1,t2,len,xnum,xint,number_1,flag_number[5]={60,60,60,60,60},jp=0,ijk=0;
   double xprob,prob;
   char sign,number_1_char[2];
   node_2=head;
   if(head==NULL)
   {
	   printf("The linked list not built successfully");
   }
   else
   {
	   for(t2=0;t2<TIMEPOINT;t2++)
		   {
		     while(t2!= node_2->time_point)
			 node_2=(LTL_Node*)node_2->nextptr;

			 for(t1=0;t1<num_prop;t1++)
			 {
               len=strlen(proposition[t1]);

				   for(ijk=0;ijk<5;ijk++)
				   {
					   flag_number[ijk]=60;
				   }
				   if(proposition[t1][0]=='P')
				   {
					      jp=0;
						  itoa(t1+1,number_1_char, 10);
						  for(number_1=0;number_1<MAX;number_1++)
						  {
							if(strcmp(node_2->LTL_Node_Element[number_1].Formula,number_1_char)==0)
							{
								flag_number[jp]=number_1;
								jp++;
							}
						  }
                       
					   xnum=atoi(&proposition[t1][2]);
					   xint=atoi(&proposition[t1][5]);
					   //xint=atoi(&proposition[t1][6]);
					   sign=proposition[t1][8];
					   //sign=proposition[t1][9];

					   xprob=atof(&proposition[t1][9]);
					   //xprob=atof(&proposition[t1][10]);
					   prob=brown[xnum][t2][xint];
					   if(sign=='<')
                           {
                                if(prob<xprob)
                                                 {
                                                      for(jp=0;jp<4;jp++)
													  {
													  if(flag_number[jp]!=60)
													  {
													  node_2->LTL_Node_Element[flag_number[jp]].Truth_flag=1;
                                                      node_2->LTL_Node_Element[flag_number[jp]].evaluated=1;
													  }
													  }
                                                           
                                                 }
								else if(prob>=xprob)
                                                 {
													  for(jp=0;jp<4;jp++)
													  { 
													  if(flag_number[jp]!=60)
													  {
													  node_2->LTL_Node_Element[flag_number[jp]].Truth_flag=0;
                                                      node_2->LTL_Node_Element[flag_number[jp]].evaluated=1;
													  }
													  }
                                                           
                                                 }
                             }
                        if(sign=='>')
                              {
                                if(prob>xprob)
                                                 {
                                                      for(jp=0;jp<4;jp++)
													  {
													  if(flag_number[jp]!=60)
													  {
													  node_2->LTL_Node_Element[flag_number[jp]].Truth_flag=1;
                                                      node_2->LTL_Node_Element[flag_number[jp]].evaluated=1;
													  }
													  }
                                                           
                                                 }
								else if(prob<=xprob)
                                                 {
                                                      for(jp=0;jp<4;jp++)
													  {
													  if(flag_number[jp]!=60)
													  {
													  node_2->LTL_Node_Element[flag_number[jp]].Truth_flag=0;
                                                      node_2->LTL_Node_Element[flag_number[jp]].evaluated=1;
													  }
													  }
                                                           
                                                 }
                             }
                          if(sign=='=')
                             {
                                 if(prob==xprob)
                                         {
                                            for(jp=0;jp<4;jp++)
											{   
											if(flag_number[jp]!=60)
											{
											node_2->LTL_Node_Element[flag_number[jp]].Truth_flag=1;
                                            node_2->LTL_Node_Element[flag_number[jp]].evaluated=1;
											}
											}
                                          }
								 else
								 {
									for(jp=0;jp<4;jp++)
									{
									if(flag_number[jp]!=60)
									{
									 node_2->LTL_Node_Element[flag_number[jp]].Truth_flag=0;
                                     node_2->LTL_Node_Element[flag_number[jp]].evaluated=1;
									}
									}
								 }
						    }
				   }
			 }
	   }
   }
                         

}

void update_node_remaining()
{
    int num_counter,new_num_counter,found,halftrue=0,Halftrue_2=0,t2, Formula_len,x,and_flag,and_count=0,or_count=0,Until_flag=0,i=0,openbrac_c=0,Equivalence_flag=0,Implies_flag=0;
	int or_flag=0,ijk=0,gone_in=0,len_count_2=0,global_formula_flag=0,some_count3=0,pp=-1,some_count4=0,flag_count=0,and_count2=0,len_counter_2[250]={0},int_array[10]={0},ppp,Prop_truth_flag=0,flag_number=0,evaluation_flag=0;
	char temp_var[250],test_str[10][250],test_str1[10];
	
	LTL_Node *node_2;
	node_2=tail;
	memset(temp_var,'\0',250);
	memset(test_str1,'\0',10);
	

	for(num_counter=Formula_count;num_counter>=0;num_counter--)
	{
		
		node_2=tail;
		and_count=0;
		or_count=0;
		or_flag=0;
		gone_in=0;
		len_count_2=0;
		pp=-1;
		flag_count=0;
		and_count2=0;
		Until_flag=0;
		Implies_flag=0;
		Equivalence_flag=0;
		for(i=0;i<250;i++)
		{
		len_counter_2[i]=0;
		}
		for(i=0;i<10;i++)
		{
		int_array[i]=0;
		}
		Prop_truth_flag=0;
		flag_number=0;
		i=0;
		
		for(t2=TIMEPOINT-1;t2>=0;t2--)
		   {
		     and_count=0;
			 or_count=0;
			 or_flag=0;
			 and_flag=0;
			 gone_in=0;
			 len_count_2=0;
			 pp=-1;
			 flag_count=0;
			 and_count2=0;
			 Until_flag=0;
			 Implies_flag=0;
			 Equivalence_flag=0;
			 Prop_truth_flag=0;
			 evaluation_flag=0;
			 Halftrue_2=0;
			 for(i=0;i<250;i++)
			 {
				 len_counter_2[i]=0;
			 }

			 while(t2!= node_2->time_point)
			 node_2=(LTL_Node*)node_2->previousptr;
             pp=-1;
			 for(i=0;i<10;i++)
				{
				memset(test_str[i],'\0',250);
				}
			 some_count4=0;
			 openbrac_c=0;
			 global_formula_flag=0;
			 
			 Formula_len=strlen(node_2->LTL_Node_Element[num_counter].Formula);
				 
					 	 
			 if(node_2->LTL_Node_Element[num_counter].evaluated!=1 && strlen(node_2->LTL_Node_Element[num_counter].Formula)!=0)
			 {
                for(x=0;x<Formula_len;x++)
				 {
					 if(node_2->LTL_Node_Element[num_counter].Formula[x]=='&')
					 {
						 //and_flag=1;
						 and_count=and_count++;
					 }

					 else if(node_2->LTL_Node_Element[num_counter].Formula[x]=='|')
					 {
						 //or_flag=1;
						 and_count=and_count++;
					 }

					 else if(node_2->LTL_Node_Element[num_counter].Formula[x]=='U')
					 {
						 //Until_flag=1;
						 and_count=and_count++;
					 }

					 else if(node_2->LTL_Node_Element[num_counter].Formula[x]=='-' && node_2->LTL_Node_Element[num_counter].Formula[x-1]=='<' && node_2->LTL_Node_Element[num_counter].Formula[x+1]=='>')
					 {
						//Equivalence_flag=1;
						and_count=and_count++;
					 }
					 else if(node_2->LTL_Node_Element[num_counter].Formula[x]=='-' && node_2->LTL_Node_Element[num_counter].Formula[x-1]!='<' && node_2->LTL_Node_Element[num_counter].Formula[x+1]=='>')
					 {
						//Implies_flag=1;
						and_count=and_count++;
					 }

				 }
                
				i=0;
				while(i<Formula_len)
						 {
						 if(node_2->LTL_Node_Element[num_counter].Formula[i]=='(')
								{
									for(openbrac_c=i;openbrac_c<Formula_len;openbrac_c++)
									{
										len_counter_2[openbrac_c]++;
									}
								}
								else if(node_2->LTL_Node_Element[num_counter].Formula[i]==')')
								{
									for(openbrac_c=i;openbrac_c<Formula_len;openbrac_c++)
									{
										len_counter_2[openbrac_c]--;
									}
								}
						 i++;
						 }	
				for(some_count3=0;some_count3<Formula_len;some_count3++)
				{
						if(node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='&' || node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='|' || node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='U'|| node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='-')
						{
						  if(len_counter_2[some_count3]==0)
						  {
							  global_formula_flag=1;
						  }
						}

			// added 12/10
						if(node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='&' && len_counter_2[some_count3]==0)
					 {
						 and_flag=1;
						
					 }

						else if(node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='|'&& len_counter_2[some_count3]==0)
					 {
						 or_flag=1;
						 
					 }

						else if(node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='U'&& len_counter_2[some_count3]==0)
					 {
						 Until_flag=1;
						 
					 }

						else if(node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='-' && node_2->LTL_Node_Element[num_counter].Formula[x-1]=='<' && node_2->LTL_Node_Element[num_counter].Formula[x+1]=='>'&& len_counter_2[some_count3]==0)
					 {
						 Equivalence_flag=1;
						 
					 }
						else if(node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='-' && node_2->LTL_Node_Element[num_counter].Formula[x-1]!='<' && node_2->LTL_Node_Element[num_counter].Formula[x+1]=='>'&& len_counter_2[some_count3]==0)
					 {
						 Implies_flag=1;
						 
					 }


			// added 12/10


				}
				 	 
				 
				 
				if(node_2->LTL_Node_Element[num_counter].Formula[0]=='G' && global_formula_flag!=1)
				{
                   
				   memset(temp_var,'\0',250);
				   strncpy(temp_var,&node_2->LTL_Node_Element[num_counter].Formula[2],strlen(node_2->LTL_Node_Element[num_counter].Formula)-3); 
				   found=0;
				   while(found==0)
				   {
					   for(new_num_counter=num_counter+1;new_num_counter<Formula_count;new_num_counter++)
					   {
						   if(strcmp(temp_var,node_2->LTL_Node_Element[new_num_counter].Formula)==0 && node_2->LTL_Node_Element[new_num_counter].evaluated ==1 )
						   {
							   found=1;
							   if(node_2->LTL_Node_Element[new_num_counter].Truth_flag==1)
							   {
								   halftrue=1;
							   }
							   else if(node_2->LTL_Node_Element[new_num_counter].Truth_flag==0)
							   {
                                  node_2->LTL_Node_Element[num_counter].Truth_flag=0;
								  node_2->LTL_Node_Element[num_counter].evaluated=1;
							   }

                               if(halftrue==1)
							   {
								   if(t2==TIMEPOINT-1)
								   {
                                      node_2->LTL_Node_Element[num_counter].Truth_flag=1;
								      node_2->LTL_Node_Element[num_counter].evaluated=1;
								   }
								   else
								   {
								   node_2=(LTL_Node*)node_2->nextptr;
								   if(node_2->LTL_Node_Element[num_counter].Truth_flag==1)
								   {
								      node_2=(LTL_Node*)node_2->previousptr;
                                      node_2->LTL_Node_Element[num_counter].Truth_flag=1;
								      node_2->LTL_Node_Element[num_counter].evaluated=1;
								   }
								   if(node_2->LTL_Node_Element[num_counter].Truth_flag==0)
								   {
								      node_2=(LTL_Node*)node_2->previousptr;
                                      node_2->LTL_Node_Element[num_counter].Truth_flag=0;
								      node_2->LTL_Node_Element[num_counter].evaluated=1;
								   }
								   }
							   }
						   }
					   }
				   }
				}


				if(node_2->LTL_Node_Element[num_counter].Formula[0]=='X'&& global_formula_flag!=1)
								{
								   memset(temp_var,'\0',250); 
								   strncpy(temp_var,&node_2->LTL_Node_Element[num_counter].Formula[2],strlen(node_2->LTL_Node_Element[num_counter].Formula)-3); 
								   
								   found=0;
								   while(found==0)
								   {
									   for(new_num_counter=num_counter+1;new_num_counter<Formula_count;new_num_counter++)
									   {
										   if(strcmp(temp_var,node_2->LTL_Node_Element[new_num_counter].Formula)==0 && node_2->LTL_Node_Element[new_num_counter].evaluated ==1 )
										   {
											   found=1;
											     if(t2==TIMEPOINT-1)
												   {
													  node_2->LTL_Node_Element[num_counter].Truth_flag=0;
													  node_2->LTL_Node_Element[num_counter].evaluated=1;
												   }
												 else
												   {
												   node_2=(LTL_Node*)node_2->nextptr;
												   if(node_2->LTL_Node_Element[num_counter].Truth_flag==1)
												   {
													  node_2=(LTL_Node*)node_2->previousptr;
													  node_2->LTL_Node_Element[num_counter].Truth_flag=1;
													  node_2->LTL_Node_Element[num_counter].evaluated=1;
												   }
												     if(node_2->LTL_Node_Element[num_counter].Truth_flag==0)
												   {
													  node_2=(LTL_Node*)node_2->previousptr;
													  node_2->LTL_Node_Element[num_counter].Truth_flag=0;
													  node_2->LTL_Node_Element[num_counter].evaluated=1;
												   }
												   }
											   }
										   }
									   }
				                }
				if(node_2->LTL_Node_Element[num_counter].Formula[0]=='F'&& global_formula_flag!=1)
				{
                   memset(temp_var,'\0',250); 
                   strncpy(temp_var,&node_2->LTL_Node_Element[num_counter].Formula[2],strlen(node_2->LTL_Node_Element[num_counter].Formula)-3); 
				   
				   found=0;
				   while(found==0)
				   {
					   for(new_num_counter=num_counter+1;new_num_counter<Formula_count;new_num_counter++)
					   {
						   if(strcmp(temp_var,node_2->LTL_Node_Element[new_num_counter].Formula)==0 && node_2->LTL_Node_Element[new_num_counter].evaluated ==1 )
						   {
							   found=1;
							   if(node_2->LTL_Node_Element[new_num_counter].Truth_flag==1)
							   {
								      node_2->LTL_Node_Element[num_counter].Truth_flag=1;
								      node_2->LTL_Node_Element[num_counter].evaluated=1;
							   }
							   else if(node_2->LTL_Node_Element[new_num_counter].Truth_flag==0)
							   {
                                 if(t2==TIMEPOINT-1)
										{
													  node_2->LTL_Node_Element[num_counter].Truth_flag=0;
													  node_2->LTL_Node_Element[num_counter].evaluated=1;
										}
								 else
								 {
								   node_2=(LTL_Node*)node_2->nextptr;
								   if(node_2->LTL_Node_Element[num_counter].Truth_flag==1)
								   {
								      node_2=(LTL_Node*)node_2->previousptr;
                                      node_2->LTL_Node_Element[num_counter].Truth_flag=1;
								      node_2->LTL_Node_Element[num_counter].evaluated=1;
								   }
								  else if(node_2->LTL_Node_Element[num_counter].Truth_flag==0)
								   {
								      node_2=(LTL_Node*)node_2->previousptr;
                                      node_2->LTL_Node_Element[num_counter].Truth_flag=0;
								      node_2->LTL_Node_Element[num_counter].evaluated=1;
								   }
								 }
							   }

                              }
						   }
					   }
				   }


				if(node_2->LTL_Node_Element[num_counter].Formula[0]=='!'&& global_formula_flag!=1)
				{
                   memset(temp_var,'\0',250); 
					   if(node_2->LTL_Node_Element[num_counter].Formula[1]=='(')
					   {
						   strncpy(temp_var,&node_2->LTL_Node_Element[num_counter].Formula[2],strlen(node_2->LTL_Node_Element[num_counter].Formula)-3); 
					   }
					   else
					   {
						   //if(strlen(&node_2->LTL_Node_Element[num_counter].Formula)==2)
						   //{
						   strncpy(temp_var,&node_2->LTL_Node_Element[num_counter].Formula[1],strlen(node_2->LTL_Node_Element[num_counter].Formula)-1); 
						   //}
					   }
				   found=0;
				   while(found==0)
				   {
					   for(new_num_counter=num_counter+1;new_num_counter<Formula_count;new_num_counter++)
					   {
						   if(strcmp(temp_var,node_2->LTL_Node_Element[new_num_counter].Formula)==0 && node_2->LTL_Node_Element[new_num_counter].evaluated ==1 )
						   {
							   found=1;
							   if(node_2->LTL_Node_Element[new_num_counter].Truth_flag==1)
							   {
								      node_2->LTL_Node_Element[num_counter].Truth_flag=0;
								      node_2->LTL_Node_Element[num_counter].evaluated=1;
							   }
							   else if(node_2->LTL_Node_Element[new_num_counter].Truth_flag==0)
							   {
								   node_2->LTL_Node_Element[num_counter].Truth_flag=1;
								   node_2->LTL_Node_Element[num_counter].evaluated=1;
                                 
							   }

                              }
						   }
					   }
				   }
               
                 
                 i=0;
				 if(global_formula_flag==1)
				 {
									 
				 for(some_count3=0;some_count3<Formula_len;some_count3++)
						{
						if(node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='&' || node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='|'||node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='U')
						{
						   and_count2=and_count2+1;
						   if(len_counter_2[some_count3]==0)
						   {
							   pp=pp+1;
							   strncpy(test_str[pp],&node_2->LTL_Node_Element[num_counter].Formula[some_count4],some_count3-some_count4); 
							   some_count4= some_count3+1;
							   flag_count=1;
						  	}
						   if(and_count2==and_count && flag_count==1)
							{
								pp=pp+1;
								strncpy(test_str[pp],&node_2->LTL_Node_Element[num_counter].Formula[some_count4],Formula_len-some_count4);
							}
						  
					   }
                       if(node_2->LTL_Node_Element[num_counter].Formula[some_count3]=='-')
					   {
                          and_count2=and_count2+1;
								   if(node_2->LTL_Node_Element[num_counter].Formula[some_count3-1]!='<' && node_2->LTL_Node_Element[num_counter].Formula[some_count3+1]=='>')
								  {
								   if(len_counter_2[some_count3]==0)
									   {
										   pp=pp+1;
										   strncpy(test_str[pp],&node_2->LTL_Node_Element[num_counter].Formula[some_count4],some_count3-some_count4); 
										   some_count4= some_count3+2;
										   flag_count=1;
						  				}
								   if(and_count2==and_count && flag_count==1)
										{
											 pp=pp+1;
											 strncpy(test_str[pp],&node_2->LTL_Node_Element[num_counter].Formula[some_count4],Formula_len-some_count4);
										}
								  }//end of in_str[k][some_count3-1]!='<' && in_str[k][some_count3-1]!='>'

								  else if(node_2->LTL_Node_Element[num_counter].Formula[some_count3-1]=='<' && node_2->LTL_Node_Element[num_counter].Formula[some_count3+1]=='>')
								  {
								   if(len_counter_2[some_count3]==0)
									   {
										   pp=pp+1;
										   strncpy(test_str[pp],&node_2->LTL_Node_Element[num_counter].Formula[some_count4],some_count3-some_count4-1); 
										   some_count4= some_count3+2;
										   flag_count=1;
						  				}
								   if(and_count2==and_count && flag_count==1)
										{
											 pp=pp+1;
											 strncpy(test_str[pp],&node_2->LTL_Node_Element[num_counter].Formula[some_count4],Formula_len-some_count4);
										}
								  }//end of in_str[k][some_count3-1]=='<' && in_str[k][some_count3-1]!='>'

						 }

        	       }//End of for(some_count3=0;some_count3<Formula_len;some_count3++)
				 pp=0;
				 for(ijk=0;ijk<10;ijk++)
				 {
					 int_array[ijk]=0;
				 }
				 
			
				 while(strlen(test_str[pp])>0)
				 {
					 if(test_str[pp][0]=='(' && test_str[pp][strlen(test_str[pp])-1]==')')
					 {
						    strncpy(test_str1,&test_str[pp][1],strlen(test_str[pp])-2); 
							memset(test_str[pp],'\0',10);
							strcpy(test_str[pp],test_str1);
							memset(test_str1,'\0',10);
					 }


					 found=0;   
					 for(new_num_counter=num_counter+1;new_num_counter<Formula_count;new_num_counter++)
					       {
						    if(found==0)
							{
							if(strcmp(test_str[pp],node_2->LTL_Node_Element[new_num_counter].Formula)==0 && node_2->LTL_Node_Element[new_num_counter].evaluated ==1 )
						    {
						        found=1;
								if(node_2->LTL_Node_Element[new_num_counter].Truth_flag==1)
							        { 
										int_array[pp]=1;
																			
								     }
								  else if(node_2->LTL_Node_Element[new_num_counter].Truth_flag==0 && gone_in!=1)
									  {
										  int_array[pp]=0;
									  }
							  }
							}
						 }
					
					 pp++;
				 }// end of while(strlen(test_str[pp])>0)

				 // if the symbol is | or &
                 if((and_flag==1 ||or_flag==1)&& evaluation_flag==0)
				 {
				 gone_in=0;
				 evaluation_flag=1;

				 if(and_flag==1)
				 {
					 for(ppp=0;ppp<pp;ppp++)
					 {
						 if(int_array[ppp]==1)
						 {
							 Prop_truth_flag=1;
						 }
						 else if(int_array[ppp]==0)
						 {
							 Prop_truth_flag=0;
							 gone_in=1;
						 }
					 }
				 }
				 if(or_flag==1)
				 {
					 for(ppp=0;ppp<pp;ppp++)
					 {
						 if(int_array[ppp]==1)
						 {
							 Prop_truth_flag=1;
						 }
						 
					 }
				 }
				                  
				 if (Prop_truth_flag==1 && gone_in!=1)
				 {
					node_2->LTL_Node_Element[num_counter].Truth_flag=1;
					node_2->LTL_Node_Element[num_counter].evaluated=1;
				 }
				 else if(Prop_truth_flag==0 || gone_in==1)
				 {
                   	node_2->LTL_Node_Element[num_counter].Truth_flag=0;
					node_2->LTL_Node_Element[num_counter].evaluated=1;
				 }
				}// end of if(and_flag==1 ||or_flag==1)

              // if the symbol is until 
			  if(Until_flag==1 && evaluation_flag==0)
			   {
                  evaluation_flag=1;
				  if(t2==TIMEPOINT-1)
				  {
                     
					  if(int_array[1]==1)
					  {
						 node_2->LTL_Node_Element[num_counter].Truth_flag=1;
						 node_2->LTL_Node_Element[num_counter].evaluated=1;
					  }
					  else if(int_array[1]==0)
					  {
                         node_2->LTL_Node_Element[num_counter].Truth_flag=0;
						 node_2->LTL_Node_Element[num_counter].evaluated=1;
					  }
				  }

				  else
				  {
                     if(int_array[1]==1)
					  {
						 node_2->LTL_Node_Element[num_counter].Truth_flag=1;
						 node_2->LTL_Node_Element[num_counter].evaluated=1;
					  }

					 else if(int_array[1]==0)
					 {
						 if(int_array[0]==1)
						 {
							 Halftrue_2=1;
						 }

						 else if(int_array[0]==0)
						 {
							 node_2->LTL_Node_Element[num_counter].Truth_flag=0;
							 node_2->LTL_Node_Element[num_counter].evaluated=1;
						 }

						 if(Halftrue_2==1)
						 {
							 node_2=(LTL_Node*)node_2->nextptr;
								   if(node_2->LTL_Node_Element[num_counter].Truth_flag==1)
								   {
								      node_2=(LTL_Node*)node_2->previousptr;
                                      node_2->LTL_Node_Element[num_counter].Truth_flag=1;
								      node_2->LTL_Node_Element[num_counter].evaluated=1;
								   }
								   if(node_2->LTL_Node_Element[num_counter].Truth_flag==0)
								   {
								      node_2=(LTL_Node*)node_2->previousptr;
                                      node_2->LTL_Node_Element[num_counter].Truth_flag=0;
								      node_2->LTL_Node_Element[num_counter].evaluated=1;
								   }
						 }//end of if(Halftrue_2==1)
					 }//end of else if(int_array[1]==0)
				  }// end of else
			   }//end  of if(Until_flag==1)

			  //if the symbol is <->

			  if(Equivalence_flag==1 && evaluation_flag==0)
			  {
				  evaluation_flag=1;
				  if(int_array[0]==0)
				  {
                     if(int_array[1]==0)
					 {
                       node_2->LTL_Node_Element[num_counter].Truth_flag=1;
					   node_2->LTL_Node_Element[num_counter].evaluated=1;
					 }

					 else if(int_array[1]==1)
					 {
                       node_2->LTL_Node_Element[num_counter].Truth_flag=0;
					   node_2->LTL_Node_Element[num_counter].evaluated=1;
					 }
				  }

				  else if(int_array[0]==1)
				  {
                     if(int_array[1]==0)
					 {
                       node_2->LTL_Node_Element[num_counter].Truth_flag=0;
					   node_2->LTL_Node_Element[num_counter].evaluated=1;
					 }

					 else if(int_array[1]==1)
					 {
                       node_2->LTL_Node_Element[num_counter].Truth_flag=1;
					   node_2->LTL_Node_Element[num_counter].evaluated=1;
					 }
				  }
			  }//End of if(Equivalence_flag==1)

			  //if the symbol is ->

			  if(Implies_flag==1 && evaluation_flag==0)
			  {
				  evaluation_flag=1;
				  if(int_array[0]==0)
				  {
                     if(int_array[1]==0)
					 {
                       node_2->LTL_Node_Element[num_counter].Truth_flag=1;
					   node_2->LTL_Node_Element[num_counter].evaluated=1;
					 }

					 else if(int_array[1]==1)
					 {
                       node_2->LTL_Node_Element[num_counter].Truth_flag=1;
					   node_2->LTL_Node_Element[num_counter].evaluated=1;
					 }
				  }

				  else if(int_array[0]==1)
				  {
                     if(int_array[1]==0)
					 {
                       node_2->LTL_Node_Element[num_counter].Truth_flag=0;
					   node_2->LTL_Node_Element[num_counter].evaluated=1;
					 }

					 else if(int_array[1]==1)
					 {
                       node_2->LTL_Node_Element[num_counter].Truth_flag=1;
					   node_2->LTL_Node_Element[num_counter].evaluated=1;
					 }
				  }
			  }//End of if(Implies_flag==1)
     
			}// if(global_formula_flag==1)

				}
								
			 }
		}
	}

//function that defines that sorting should be based on length of the string
int cstring_cmp(const void *a, const void *b) 
{ 
	return ((strlen(b)-strlen(a)));
} 


// makes a new string array with no duplicates.
char** unique_copy(char array1[][250], int number) //sorted array
{
	int k = 0, i=0;
	char **array2;
	int unique_e=0;

	for (i = 0; i < number; i++) 
	{
		for (k=i+1;k<number;k++)
		{
			if (strcmp(array1[i],array1[k])==0) 
			{ 
				if (strcmp(array1[i],"")!=0)
				{
					memset(array1[k],'\0',250);
				}
			}
		}

	}

	for (i = 0; i < number; i++) 
	{
		if (strcmp(array1[i],"")!=0)
		{
			unique_e++;
		}
	}


	array2 = (char**)malloc(unique_e*sizeof(char*));
	for (i = 0; i <unique_e; i++)
		array2[i] = (char*) malloc(250*sizeof(char));

	for (i=0;i<unique_e;i++)
	{
		memset(array2[i],'\0',250);
	}

	k=0;
	for (i = 0; i < number; i++) 
	{
		if (strcmp(array1[i],"")!=0)
		{
			strcpy(array2[k],array1[i]);
			k++;

		}
	}
	Formula_count=unique_e;
	return array2;
}
	