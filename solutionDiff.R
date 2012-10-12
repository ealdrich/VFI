  # Determine which methods to compare from command line args
  Methods = commandArgs(TRUE);
  Method1 = Methods[1];
  Method2 = Methods[2];

print(Method1)

  # Names of value and policy function files
  vDat = '/valueFunc.dat';
  pDat = '/policyFunc.dat';

  # Import value functions and compute difference
  V1 = scan(paste(Method1, vDat, sep=''), skip=2);
  V2 = scan(paste(Method2, vDat, sep=''), skip=2);
  VDiff = max(abs(V1 - V2));
  
  # Import value functions and compute difference
  G1 = scan(paste(Method1, pDat, sep=''), skip=2);
  G2 = scan(paste(Method2, pDat, sep=''), skip=2);
  GDiff = max(abs(G1 - G2));

  # Write out
  fileName = paste('Errors_', Method1, '_', Method2, '.dat', sep='');
  write(c(VDiff, GDiff), fileName, ncolumns=1);

