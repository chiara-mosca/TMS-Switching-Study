# Install and load required packages  
if (!require("DiagrammeR")) install.packages("DiagrammeR")  
if (!require("DiagrammeRsvg")) install.packages("DiagrammeRsvg")  
if (!require("rsvg")) install.packages("rsvg")

library(DiagrammeR)  
library(DiagrammeRsvg)  
library(rsvg)

# Create the participant flow diagram  
flow <- grViz("  
digraph participant_flow {

  graph [layout = dot, rankdir = TB, fontsize = 12, fontname = 'Arial']  
    
  node [shape = box, style = filled, fillcolor = '#E8E8E8',   
        fontname = 'Arial', fontsize = 10, width = 4, margin = 0.2]

  # Level 1: REDCap export  
  A [label = 'REDCap Export\\nN = 146 patients\\n(1,211 rows × 132 columns)']

  # Level 2: Remove test row  
  B [label = 'After removing test/template row\\nN = 145 patients']  
    
  exc1 [label = 'Excluded: n = 1\\n• Test/template row (TMS001)',   
        fillcolor = '#FFD0D0', shape = box]

  # Level 3: Remove no PHQ-9 data  
  C [label = 'After initial exclusions\\nN = 130 patients']  
    
  exc2 [label = 'Excluded: n = 15\\n• No usable PHQ-9 data\\n(no forms in folder, rTMS study\\nparticipants with blank course 1,\\ndiscontinued courses with\\nno assessments)',   
        fillcolor = '#FFD0D0', shape = box]

  # Level 4: Filter to course 1 + clean  
  D [label = 'After filtering to Course 1 only\\n& removing invalid/missing PHQ-9 data\\nN = 113 patients']  
    
  exc3 [label = 'Excluded: n = 17\\n• Multiple course filtering\\n• Invalid PHQ-9 data\\n• Missing PHQ-9 data',   
        fillcolor = '#FFD0D0', shape = box]

  # Level 5: Protocol restrictions  
  E [label = 'Final Analytic Sample\\nN = 85 patients',   
     fillcolor = '#C8E6C9']  
    
  exc4 [label = 'Excluded: n = 28\\n• 26 used non-iTBS protocols\\n   (10 Hz, 10 Hz/1 Hz, or mixed)\\n• 1 switched bilateral → iTBS\\n   (wrong direction; TMS093)\\n• 1 alternated between protocols\\n   (no clean switch; TMS084)',   
        fillcolor = '#FFD0D0', shape = box]

  # Level 6: Group allocation  
  F [label = 'iTBS-Only Group\\nn = 41', fillcolor = '#BBDEFB']  
  G [label = 'Switcher Group\\n(iTBS → Bilateral)\\nn = 44', fillcolor = '#BBDEFB']

  # Edges  
  A -> B  
  A -> exc1 [style = dashed]  
  B -> C  
  B -> exc2 [style = dashed]  
  C -> D  
  C -> exc3 [style = dashed]  
  D -> E  
  D -> exc4 [style = dashed]  
  E -> F  
  E -> G

  # Alignment  
  {rank = same; A}  
  {rank = same; B; exc1}  
  {rank = same; C; exc2}  
  {rank = same; D; exc3}  
  {rank = same; E; exc4}  
  {rank = same; F; G}  
}  
")

# Display the diagram  
flow

# Save as PNG  
flow_svg <- export_svg(flow)  
rsvg_png(charToRaw(flow_svg), file = "participant_flow_diagram.png", width = 1200, height = 1600)

# Save as PDF  
rsvg_pdf(charToRaw(flow_svg), file = "participant_flow_diagram.pdf", width = 8, height = 11)

cat("Flow diagram saved as participant_flow_diagram.png and participant_flow_diagram.pdf\n")  