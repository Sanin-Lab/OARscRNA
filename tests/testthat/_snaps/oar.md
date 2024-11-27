# oar works with defaults

    Code
      waldo::compare(oar_output, check_oar_meta)
    Output
      `old` is length 13
      `new` is length 38
      
          names(old)     | names(new)                    
      [5] "condition"    | "condition"    [5]            
      [6] "nCount_SCT"   | "nCount_SCT"   [6]            
      [7] "nFeature_SCT" | "nFeature_SCT" [7]            
                         - "res_0"        [8]            
                         - "res_0.05"     [9]            
                         - "res_0.1"      [10]           
                         - "res_0.15"     [11]           
                         - "res_0.2"      [12]           
                         - "res_0.25"     [13]           
                         - "res_0.3"      [14]           
      ... ...              ...            and 21 more ...
      
      `old$res_0` is absent
      `new$res_0` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.05` is absent
      `new$res_0.05` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.1` is absent
      `new$res_0.1` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.15` is absent
      `new$res_0.15` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.2` is absent
      `new$res_0.2` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.25` is absent
      `new$res_0.25` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.3` is absent
      `new$res_0.3` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.35` is absent
      `new$res_0.35` is an S3 object of class <factor>, an integer vector
      
      And 17 more differences ...

# oar works with no mismatch

    Code
      waldo::compare(oar_output, check_oar_no_mis_meta)
    Output
      `old` is length 13
      `new` is length 38
      
          names(old)     | names(new)                    
      [5] "condition"    | "condition"    [5]            
      [6] "nCount_SCT"   | "nCount_SCT"   [6]            
      [7] "nFeature_SCT" | "nFeature_SCT" [7]            
                         - "res_0"        [8]            
                         - "res_0.05"     [9]            
                         - "res_0.1"      [10]           
                         - "res_0.15"     [11]           
                         - "res_0.2"      [12]           
                         - "res_0.25"     [13]           
                         - "res_0.3"      [14]           
      ... ...              ...            and 21 more ...
      
      `old$res_0` is absent
      `new$res_0` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.05` is absent
      `new$res_0.05` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.1` is absent
      `new$res_0.1` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.15` is absent
      `new$res_0.15` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.2` is absent
      `new$res_0.2` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.25` is absent
      `new$res_0.25` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.3` is absent
      `new$res_0.3` is an S3 object of class <factor>, an integer vector
      
      `old$res_0.35` is absent
      `new$res_0.35` is an S3 object of class <factor>, an integer vector
      
      And 17 more differences ...

