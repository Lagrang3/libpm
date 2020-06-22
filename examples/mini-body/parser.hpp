#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

class CmdParser
{
   private:
    std::string executable;
    std::set<std::string> opt_noarg;
    std::map<std::string, std::string> opt_warg;

   public:
    CmdParser(int argc, char** argv) : executable{argv[0]}
    {
        for (int i = 1; i < argc;)
        {
            std::string opt{argv[i++]};

            if (opt.size() == 2 and opt[0] == '-')
            {
                // short options
                opt = opt.substr(1);
            }
            else if (opt.size() > 2 and opt.substr(0, 2) == "--")
            {
                // long options
                opt = opt.substr(2);
            }
            else
            {
                throw std::runtime_error("Invalid option or command: " + opt);
            }

            if (i >= argc)
            {
                continue;
                opt_noarg.insert(opt);
            }
            std::string next{argv[i]};
            if (next[0] == '-')
            {
                opt_noarg.insert(opt);
            }
            else
            {
                opt_warg[opt] = next;
                ++i;
            }
        }
    }

    auto& operator[](const std::string& key) const { return opt_warg.at(key); }
    bool find(const std::string& key) const
    {
        return opt_noarg.find(key) != opt_noarg.end();
    }
};

class ParamfileParser
{
   private:
    std::map<std::string, std::string> options;

   public:
    ParamfileParser(const std::string& filename)
    {
        std::ifstream input(filename);
        for (std::string line; std::getline(input, line);)
        {
            line.erase(std::remove(line.begin(), line.end(), ' '),
                       line.end());  // remove spaces
            line.erase(std::remove(line.begin(), line.end(), '\t'),
                       line.end());  // remove tabs
            line.erase(std::find(line.begin(), line.end(), '#'),
                       line.end());  // remove comments

            if (line == "")
                continue;
            auto middle =
                std::find(line.begin(), line.end(), '=') - line.begin();
            std::string left = line.substr(0, middle),
                        right = line.substr(middle + 1);

            options[left] = right;
        }
    }
    const auto& operator[](const std::string& key) const
    {
        return options.at(key);
    }
};
