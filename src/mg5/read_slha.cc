#include "tthAnalysis/tthMEM/interface/mg5/read_slha.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line_ext()

#include <algorithm>
#include <fstream>

using namespace std;

void SLHABlock::set_entry(vector<int> indices, double value)
{
  if (_entries.size() == 0)
    _indices = indices.size();
  else if(indices.size() != _indices)
    throw_line_ext("runtime error", TTHEXCEPTION_ERR_CODE_SLHA)
      << "Wrong number of indices in set_entry";

  _entries[indices] = value;
}

double SLHABlock::get_entry(vector<int> indices, double def_val)
{
  if (_entries.find(indices) == _entries.end()){
    LOGWARN << "Warning: No such entry in " << _name << ", using default value "
         << def_val;
    return def_val;
  }
  return _entries[indices];
}

void SLHAReader::read_slha_file(string file_name)
{
  ifstream param_card;
  param_card.open(file_name.c_str(), ifstream::in);
  if(!param_card.good())
    throw_line_ext("input error", TTHEXCEPTION_ERR_CODE_FILE_NOT_FOUND)
      << "Error while opening param card";
  LOGINFO << "Opened slha file " << file_name << " for reading";
  char buf[200];
  string line;
  string block("");

  while(param_card.good()){
    param_card.getline(buf, 200);
    line = buf;
    // Change to lowercase
    transform(line.begin(), line.end(), line.begin(), (int(*)(int)) tolower);
    if(line != "" && line[0] != '#'){
      if(block != ""){
	// Look for double index blocks
	double dindex1, dindex2;
	double value;
	stringstream linestr2(line);
	if (linestr2 >> dindex1 >> dindex2 >> value &&
	    dindex1 == int(dindex1) and dindex2 == int(dindex2))
	{
	  vector<int> indices;
	  indices.push_back(int(dindex1));
	  indices.push_back(int(dindex2));
	  set_block_entry(block, indices, value);
	  // Done with this line, read next
	  continue;
	}
	stringstream linestr1(line);
	// Look for single index blocks
	if(linestr1 >> dindex1 >> value && dindex1 == int(dindex1)) 
	{
	  vector<int> indices;
	  indices.push_back(int(dindex1));
	  set_block_entry(block, indices, value);
	  // Done with this line, read next
	  continue;
	}
      }
      // Look for block
      if(line.find("block ") != line.npos){
	line = line.substr(6);
	// Get rid of spaces between block and block name
	while (line[0] == ' ')
	  line = line.substr(1);
	// Now find end of block name
	size_t space_pos = line.find(' ');
	if(space_pos != line.npos)
	  line = line.substr(0, space_pos);
	block = line;
	continue;
      }
      // Look for decay
      if(line.find("decay ") == 0){
	line = line.substr(6);
	block = "";
	stringstream linestr(line);
	int pdg_code;
	double value;
	if(linestr >> pdg_code >> value)
	  set_block_entry("decay", pdg_code, value);
	else
	  LOGWARN << "Warning: Wrong format for decay block " << line;
	continue;
      }
    }
  }

  if (_blocks.size() == 0)
    throw_line_ext("runtime error", TTHEXCEPTION_ERR_CODE_SLHA)
      << "No information read from SLHA card";

  param_card.close();
}

double SLHAReader::get_block_entry(string block_name, vector<int> indices, 
				   double def_val)
{
  if (_blocks.find(block_name) == _blocks.end()){
    LOGWARN << "No such block " << block_name << ", using default value "
         << def_val;
    return def_val;
  }
  return _blocks[block_name].get_entry(indices);  
}

double SLHAReader::get_block_entry(string block_name, int index, 
				   double def_val)
{
  vector<int> indices;
  indices.push_back(index);
  return get_block_entry(block_name, indices, def_val);
}


void SLHAReader::set_block_entry(string block_name, vector<int> indices, 
				   double value)
{
  if (_blocks.find(block_name) == _blocks.end()){
    SLHABlock block(block_name);
    _blocks[block_name] = block;
  }
  _blocks[block_name].set_entry(indices, value);  
  /* cout << "Set block " << block_name << " entry ";
     for (int i=0;i < indices.size();i++) 
     cout << indices[i] << ' ';
     cout << "to " << _blocks[block_name].get_entry(indices) << '\n';*/
}

void SLHAReader::set_block_entry(string block_name, int index, 
				   double value)
{
  vector<int> indices;
  indices.push_back(index);
  set_block_entry(block_name, indices, value);
}
