MODULE Main;

IMPORT LR3, Simulator, KinematicConstraint;
FROM Simulator IMPORT AddConstraint;

BEGIN
  Simulator.Init();
  Simulator.Go();
END Main.
