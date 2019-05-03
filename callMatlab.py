## JAVA CODE
# public static void callMatlab(ArrayList<NmrStructure> structures) throws 
#   	MatlabConnectionException, MatlabInvocationException {

# 	  MatlabProxyFactoryOptions options = new MatlabProxyFactoryOptions.Builder().setUsePreviouslyControlledSession(true).build();
# 	  MatlabProxyFactory factory = new MatlabProxyFactory(options);
# 	  MatlabProxy proxy = factory.getProxy();

# 	  proxy.eval("cd ..");
# 	  proxy.eval("cd NMR_java");
# 	  proxy.eval("cd matlab");
	  
# //	  proxy.eval("addpath(genpath('/spinach/'));");
# 	  // this is actually create the matlab script for the purpose
	  
# 	  for (NmrStructure nmr_str : structures) {
# 		  System.out.println(nmr_str.hmdb_id);
# 		  String isotopes = "sys.isotopes = {";
# 		  String ppms = "inter.zeeman.scalar = {";
		  
# 		  // Empty any existing coupling constants from previous runs
# 		  proxy.eval("inter.coupling.scalar = []");
		  
# 		  for (int i = 0; i < nmr_str.hydrogen_positions.size(); i++) {
# 			  if (i == nmr_str.hydrogen_positions.size() - 1) {
# 				  isotopes = isotopes + "'1H'}";
# 				  ppms = ppms + String.valueOf(nmr_str.chemical_shifts.get(i)) + "}";
# 				  String coupling = "inter.coupling.scalar{" + String.valueOf(i+1) + "," +
# 			      String.valueOf(i+1) + "} = 0" ;
# 				  System.out.println("coupling: "+ coupling);
# 				  proxy.eval(coupling);
# 			  }
# 			  else {
# 				  isotopes = isotopes + "'1H', ";
# 				  ppms = ppms + String.valueOf(nmr_str.chemical_shifts.get(i)) + ", ";
# 			  }
			  
# 		 for (int j = i + 1; j < nmr_str.hydrogen_positions.size(); j++) {
# 			 String coupling = "inter.coupling.scalar{" + String.valueOf(i+1) + "," + String.valueOf(j+1) + "} = 1" ;
# 			 System.out.println("coupling: "+coupling);
# 			 proxy.eval(coupling);
#     		}
# 		}
#     System.out.println(ppms); //	inter.zeeman.scalar = {3.07, 3.16, 3.68, 3.96, 7.0, 7.67}
#   	proxy.eval(isotopes);
#   	proxy.eval(ppms);

#   	proxy.feval("create_nmr1H_plot");
# //  
#   	System.out.println(proxy.toString());
# 	}

# 	proxy.disconnect();
# }


import matlab.engine
import os, sys

# test command: 

# call matlab script and pass argument
# Get pi from the MATLAB workspace and copy it to a Python variable.

# import matlab.engine
# eng = matlab.engine.start_matlab()
# eng.eval('a = pi;',nargout=0)
# mpi = eng.workspace['a']
# print(mpi)
# 3.14159265359
# params must include the shift for hydrogen and the hydrogen position
def callMatlab(shift,hydrogen_positions):
	
	eng = matlab.engine.start_matlab()
	eng.eval('cd matlab',nargout=0)

	isotopes = "sys.isotopes = {"
	ppms = "inter.zeeman.scalar = {"
	eng.eval("inter.coupling.scalar = []",nargout=0)		# Empty any existing coupling constants from previous runs
	for i in range(0,len(hydrogen_positions)):
		if i == (len(hydrogen_positions)-1): 
			isotopes = isotopes + "'1H'}"
			shift_value = str(shift[i])
			ppms = ppms + shift_value + "}"
			coupling = "inter.coupling.scalar{" + str(i+1) + "," + str(i+1) + "} = 0" 
			# print(coupling) # checked
			eng.eval(coupling,nargout=0)
		else:
			isotopes = isotopes + "'1H', "
			ppms = ppms + str(shift[i]) + ", "

		for j in range(i+1,len(hydrogen_positions)):
			coupling = "inter.coupling.scalar{" + str(i+1) + "," + str(j+1) + "} = 1" 
			# print(coupling) # checked
			eng.eval(coupling,nargout=0)

	# print(isotopes) # checked
	# print(ppms)     # checked
	eng.eval(isotopes,nargout=0)
	eng.eval(ppms,nargout=0)
	eng.feval("create_nmr1H_plot",nargout=0)
	eng.quit()








def main():
	shift = [3.07, 3.16, 3.68, 3.96, 7.0, 7.67]
	hydrogen_positions = [0,1,2,3,4,5]
	callMatlab(shift,hydrogen_positions)






if __name__ == '__main__':
	main()


















