    # loads in the visual c++ 2019 compiler chain
    
    cmd.exe /c "call `"C:\Devel\Lang\MSVC-2019\VC\Auxiliary\Build\vcvars64.bat`" && set > %temp%\vcvars.txt"
    
    Get-Content "$env:temp\vcvars.txt" | Foreach-Object {
      if ($_ -match "^(.*?)=(.*)$") {
        Set-Content "env:\$($matches[1])" $matches[2]
      }
    }