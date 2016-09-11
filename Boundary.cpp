Walls::Walls(int nozzleLength, const int d[2])
{
  nbNodes = 2 * d[0] + (d[1] -nozzleLength-1);
  nozzleStart = floor((d[1] - (nozzleLength + 1))*0.5);
  nodes = new int*[nbNodes];
  for (int i=0;i<nbNodes;i++)
    {nodes[i] = new int[2];}

  int cc = 0;
  for (int x=0;x<d[0];x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = d[1];
      cc += 1;
    }
  for (int x=0;x<d[0];x++)
    {
      nodes[cc][0] = x; nodes[cc][1] = 0;
      cc += 1;
    }
  for (int y=0;y<nozzleStart;y++)
    {
      nodes[cc][0] = 0;  nodes[cc][1] = y;
      cc += 1;
    }
  for (int y=nozzleStart+nozzleLength;y<d[1];y++)
    {
      nodes[cc][0] = 0; nodes[cc][1] = y;
      cc += 1;
    }
}
virtual void Walls::boundaryCondition(LatticeSite **sites, LatticeSite **_sites)
{
  int x,y;

  for (int i=0;i<nbNodes;i++)
    {
      x = nodes[i][0]; y = nodes[i][1];
      for (int k=0;k<9;k++)
	{
	  _sites[x][y].f[k] = sites[x][y].f[op[k]];
	}
    }
}
