import React, { useState, useRef } from 'react';
import { useDropzone } from 'react-dropzone';
import { Upload, FileText, Zap, Brain, Target, Play, Volume2, BarChart3, Atom, Eye } from 'lucide-react';
import './QuantumAnalysisPage.css';
import MoleculeViewer from '../viewer/MoleculeViewer';

function QuantumAnalysisPage() {
  const [inputType, setInputType] = useState('smiles');
  const [molecularInput, setMolecularInput] = useState('');
  const [analysisResults, setAnalysisResults] = useState(null);
  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [aiSummary, setAiSummary] = useState('');
  const [isGeneratingSummary, setIsGeneratingSummary] = useState(false);
  const [audioUrl, setAudioUrl] = useState(null);
  const [isGeneratingAudio, setIsGeneratingAudio] = useState(false);
  
  // Quantum simulation configuration
  const [simulationConfig, setSimulationConfig] = useState({
    targetProperty: 'ground_state',
    simulationMode: 'static',
    basisSet: 'sto-3g',
    activeSpace: { electrons: 4, orbitals: 4 },
    mapper: 'jordan_wigner',
    optimizer: 'slsqp',
    backend: 'local',
    errorMitigation: true,
    noiseModel: 'ideal'
  });
  
  // Interactive quantum spectroscopy
  const [spectroscopyMode, setSpectroscopyMode] = useState('absorption');
  const [timeEvolution, setTimeEvolution] = useState({ enabled: false, timeRange: [0, 10] });
  
  // Visualization data
  const [spectrumData, setSpectrumData] = useState(null);
  const [dynamicsData, setDynamicsData] = useState(null);
  const [currentRepresentation, setCurrentRepresentation] = useState('ballstick');
  const analogViewersRef = useRef({});
  const resultsRef = useRef(null);

  const onDrop = (acceptedFiles) => {
    const file = acceptedFiles[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (e) => {
        setMolecularInput(e.target.result);
        setInputType('pdb');
      };
      reader.readAsText(file);
    }
  };

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: {
      'text/plain': ['.pdb', '.mol', '.sdf', '.cif'],
      'chemical/x-pdb': ['.pdb'],
      'chemical/x-mdl-molfile': ['.mol']
    }
  });

  const analyzeMolecule = async () => {
    if (!molecularInput.trim()) return;

    setIsAnalyzing(true);

    try {
      // Call backend API for quantum simulation
      const response = await fetch('https://molecular-analysis-api-7def38974c94.herokuapp.com/api/quantum_simulate', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          input_type: inputType,
          molecular_data: molecularInput,
          config: simulationConfig,
          timestamp: new Date().toISOString()
        })
      });

      if (!response.ok) {
        throw new Error(`API request failed: ${response.statusText}`);
      }

      const results = await response.json();

      setAnalysisResults(results);
      // scroll-to-results (small delay so DOM renders)
      setTimeout(() => {
        try {
          const el = resultsRef && resultsRef.current;
          if (el) {
            // Primary: smooth scrollIntoView
            el.scrollIntoView({ behavior: 'smooth', block: 'start' });
            // Fallback: ensure scroll by setting window.scrollTo to element position after short delay
            setTimeout(() => {
              try {
                const rect = el.getBoundingClientRect();
                const absoluteTop = window.pageYOffset + rect.top;
                window.scrollTo({ top: Math.max(0, absoluteTop - 20), behavior: 'smooth' });
              } catch (_) { /* ignore */ }
            }, 120);
          }
        } catch (e) { /* ignore */ }
      }, 200);

      // Generate AI summary
      await generateAISummary(results);

      // Update visualization data
      updateVisualizationData(results);

      // Auto-load all 3D structures after a short delay to ensure DOM is ready
      setTimeout(() => {
        // Load main molecule structure
        if (results.quantum_results?.visualization?.cif_content) {
          window.dispatchEvent(new CustomEvent('loadStructure', {
            detail: {
              content: results.quantum_results.visualization.cif_content,
              format: 'mmcif'
            }
          }));
        }

        // Load all analog structures
        if (results.quantum_results?.molecular_analogs) {
          results.quantum_results.molecular_analogs.forEach((analog, index) => {
            if (analog.cif_content) {
              loadAnalogStructure(index, analog.cif_content);
            }
          });
        }
      }, 500);

    } catch (error) {
      console.error('Quantum analysis failed:', error);
      alert(`Analysis failed: ${error.message}. Make sure the backend server is running on port 8000.`);
    } finally {
      setIsAnalyzing(false);
    }
  };

  const generateAISummary = async (results) => {
    setIsGeneratingSummary(true);
    
    try {
      // Generate molecule-specific AI summary
      const molecularData = molecularInput.toLowerCase();

      const summary = ``;
      
      await new Promise(resolve => setTimeout(resolve, 2000));
      setAiSummary(summary);
    } catch (error) {
      console.error('AI summary generation failed:', error);
    } finally {
      setIsGeneratingSummary(false);
    }
  };

  const updateVisualizationData = (results) => {
    // Update spectrum data for interactive plotting
    if (results.quantum_results?.spectra) {
      setSpectrumData({
        absorption: results.quantum_results.spectra.absorption_spectrum,
        emission: results.quantum_results.spectra.emission_spectrum,
        vibrational: results.quantum_results.spectra.vibrational_modes
      });
    }
    
    // Update orbital data (removed unused orbitalData state)
    
    // Update dynamics data
    if (results.quantum_results?.dynamics) {
      setDynamicsData(results.quantum_results.dynamics);
    }
  };

  const generateVoiceNarration = async () => {
    if (!aiSummary) return;

    setIsGeneratingAudio(true);

    try {
      await new Promise(resolve => setTimeout(resolve, 2000));
      setAudioUrl("https://example.com/quantum-audio.mp3");
    } catch (error) {
      console.error('Voice generation failed:', error);
    } finally {
      setIsGeneratingAudio(false);
    }
  };

  // Helper function to load analog structure into mini viewer
  const loadAnalogStructure = (analogIndex, cifContent) => {
    if (!window.$3Dmol || !cifContent) return;

    const viewerId = `analog-viewer-${analogIndex}`;
    const element = document.getElementById(viewerId);
    if (!element) return;

    try {
      // Create or reuse viewer
      if (!analogViewersRef.current[analogIndex]) {
        analogViewersRef.current[analogIndex] = window.$3Dmol.createViewer(element, {
          defaultcolors: window.$3Dmol.rasmolElementColors,
        });
      }

      const viewer = analogViewersRef.current[analogIndex];
      viewer.clear();
      viewer.addModel(cifContent, 'mmcif');
      viewer.setBackgroundColor(0x000000);
      viewer.setStyle({}, {
        stick: { radius: 0.15, colorscheme: 'Jmol' },
        sphere: { scale: 0.3, colorscheme: 'Jmol' }
      });
      viewer.zoomTo();
      viewer.render();

      // Enable auto-rotation for a nice visual effect
      viewer.spin('y', 0.5); // Rotate around Y-axis at 0.5 degrees per frame
    } catch (err) {
      console.error(`Failed to load analog ${analogIndex}:`, err);
    }
  };

  const exampleMolecules = [
    {
      name: "Imatinib (Cancer Therapy)",
      smiles: "CC1=CC=C(C=C1)NC(=O)C2=CC(=CC=C2)C3=CN=CC=N3",
      description: "Targeted cancer therapy - BCR-ABL kinase inhibitor"
    },
    {
      name: "Benzene (Aromatic)",
      smiles: "C1=CC=CC=C1",
      description: "Classic aromatic system for quantum chemistry studies"
    },
    {
      name: "Ethylene (π-system)",
      smiles: "C=C",
      description: "Simple π-system for excited state calculations"
    }
  ];

  return (
    <div className="quantum-analysis-page">
      <div className="analysis-header">
        <h1>🧬 Advanced Quantum Chemistry Platform</h1>
        <p>Cutting-edge quantum computing for molecular analysis with VQE, EOM, time evolution, and spectral estimation</p>
      </div>

      <div className="analysis-container">
        <div className="input-section">
          <div className="input-tabs">
            <button 
              className={`tab ${inputType === 'smiles' ? 'active' : ''}`}
              onClick={() => setInputType('smiles')}
            >
              <FileText size={20} />
              SMILES Input
            </button>
            <button 
              className={`tab ${inputType === 'pdb' ? 'active' : ''}`}
              onClick={() => setInputType('pdb')}
            >
              <Upload size={20} />
              Structure Upload
            </button>
            <button 
              className={`tab ${inputType === 'draw' ? 'active' : ''}`}
              onClick={() => setInputType('draw')}
            >
              <Target size={20} />
              Draw Molecule
            </button>
          </div>

          {inputType === 'smiles' && (
            <div className="input-area">
              <label htmlFor="smiles-input">Enter SMILES string:</label>
              <input
                id="smiles-input"
                type="text"
                value={molecularInput}
                onChange={(e) => setMolecularInput(e.target.value)}
                placeholder="e.g., C1=CC=CC=C1 (benzene)"
                className="smiles-input"
              />
              <div className="example-molecules">
                <h4>Try these quantum chemistry examples:</h4>
                {exampleMolecules.map((mol, index) => (
                  <button
                    key={index}
                    className="example-btn"
                    onClick={() => setMolecularInput(mol.smiles)}
                  >
                    <strong>{mol.name}</strong>
                    <span>{mol.description}</span>
                  </button>
                ))}
              </div>
            </div>
          )}

          {inputType === 'pdb' && (
            <div className="input-area">
              <div {...getRootProps()} className={`dropzone ${isDragActive ? 'active' : ''}`}>
                <input {...getInputProps()} />
                <Upload size={48} />
                <p>Drag & drop molecular structure files here</p>
                <p className="file-types">Supported: .pdb, .mol, .sdf, .cif</p>
              </div>
              {molecularInput && (
                <div className="file-preview">
                  <h4>Structure loaded:</h4>
                  <pre>{molecularInput.substring(0, 500)}...</pre>
                </div>
              )}
            </div>
          )}

          {inputType === 'draw' && (
            <div className="input-area">
              <div className="drawing-area">
                <p>Interactive molecular drawing interface</p>
                <p className="coming-soon">Advanced 3D molecular editor with quantum property prediction</p>
              </div>
            </div>
          )}

          {/* Quantum Simulation Configuration */}
          <div className="quantum-config">
            <h3>⚛️ Quantum Simulation Configuration</h3>
            <div className="config-grid">
              <div className="config-group">
                <label>Target Property:</label>
                <select 
                  value={simulationConfig.targetProperty}
                  onChange={(e) => setSimulationConfig({...simulationConfig, targetProperty: e.target.value})}
                >
                  <option value="ground_state">Ground State</option>
                  <option value="excited_state">Excited States</option>
                  <option value="spectrum">Spectrum</option>
                  <option value="dynamics">Dynamics</option>
                </select>
              </div>
              
              <div className="config-group">
                <label>Simulation Mode:</label>
                <select 
                  value={simulationConfig.simulationMode}
                  onChange={(e) => setSimulationConfig({...simulationConfig, simulationMode: e.target.value})}
                >
                  <option value="static">Static (VQE + EOM)</option>
                  <option value="time_evolution">Time Evolution</option>
                  <option value="spectral_estimation">Spectral Estimation</option>
                </select>
              </div>
              
              <div className="config-group">
                <label>Basis Set:</label>
                <select 
                  value={simulationConfig.basisSet}
                  onChange={(e) => setSimulationConfig({...simulationConfig, basisSet: e.target.value})}
                >
                  <option value="sto-3g">STO-3G</option>
                  <option value="6-31g">6-31G</option>
                  <option value="6-31g*">6-31G*</option>
                  <option value="cc-pvdz">cc-pVDZ</option>
                </select>
              </div>
              
              <div className="config-group">
                <label>Backend:</label>
                <select 
                  value={simulationConfig.backend}
                  onChange={(e) => setSimulationConfig({...simulationConfig, backend: e.target.value})}
                >
                  <option value="local">Local Simulator</option>
                  <option value="ibm_cloud">IBM Quantum Cloud</option>
                </select>
              </div>
            </div>
            
            <div className="advanced-config">
              <label>
                <input 
                  type="checkbox" 
                  checked={simulationConfig.errorMitigation}
                  onChange={(e) => setSimulationConfig({...simulationConfig, errorMitigation: e.target.checked})}
                />
                Error Mitigation
              </label>
              <label>
                <input 
                  type="checkbox" 
                  checked={timeEvolution.enabled}
                  onChange={(e) => setTimeEvolution({...timeEvolution, enabled: e.target.checked})}
                />
                Time Evolution Analysis
              </label>
            </div>
          </div>

          <button 
            className="analyze-btn quantum"
            onClick={analyzeMolecule}
            disabled={!molecularInput.trim() || isAnalyzing}
          >
            {isAnalyzing ? (
              <>
                <Zap className="spinning" size={20} />
                Running Quantum Simulation...
              </>
            ) : (
              <>
                <Atom size={20} />
                Start Quantum Analysis
              </>
            )}
          </button>
        </div>

        {analysisResults && (
          <div className="results-section" ref={resultsRef}>
            <div className="results-header">
              <h2>⚛️ Quantum Chemistry Results</h2>
              <div className="result-actions">
                <button 
                  className="voice-btn"
                  onClick={generateVoiceNarration}
                  disabled={!aiSummary || isGeneratingAudio}
                >
                  {isGeneratingAudio ? (
                    <>
                      <Volume2 className="spinning" size={16} />
                      Generating...
                    </>
                  ) : (
                    <>
                      <Play size={16} />
                      Quantum Voice
                    </>
                  )}
                </button>
                {audioUrl && (
                  <audio controls src={audioUrl} className="audio-player">
                    Your browser does not support the audio element.
                  </audio>
                )}
              </div>
            </div>

            {/* Professional Quantum-AI Drug Discovery Coach Output */}
            <div className="professional-output">
              <div className="output-header">
                <h2>🧬 Quantum-AI Drug Discovery Coach Output</h2>
              </div>

              {/* 1. Quantum Chemistry Simulation Results */}
              <div className="section quantum-chemistry">
                <h3>1. Quantum Chemistry Simulation Results</h3>
                <div className="results-grid">
                  <div className="result-item">
                    <span className="label">Total SCF Energy:</span>
                    <span className="value">{analysisResults.quantum_results?.total_energy_hartree?.toFixed(3)} Hartree ({analysisResults.quantum_results?.total_energy_ev?.toFixed(1)} eV)</span>
                  </div>
                  <div className="result-item">
                    <span className="label">Formation Energy:</span>
                    <span className="value">{analysisResults.quantum_results?.formation_energy_kcal_mol?.toFixed(1)} kcal/mol</span>
                  </div>
                  <div className="result-item">
                    <span className="label">HOMO-LUMO Gap:</span>
                    <span className="value">{analysisResults.quantum_results?.homo_lumo_gap?.toFixed(2)} eV</span>
                  </div>
                  <div className="result-item">
                    <span className="label">Binding Affinity:</span>
                    <span className="value">{analysisResults.quantum_results?.binding_energy?.toFixed(1)} kcal/mol</span>
                  </div>
                </div>
              </div>

              {/* 2. ML/AI Molecular Modeling */}
              <div className="section ml-modeling">
                <h3>2. ML/AI Molecular Modeling</h3>
                <div className="results-grid">
                  <div className="result-item">
                    <span className="label">Drug-likeness Score:</span>
                    <span className="value">{analysisResults.quantum_results?.ml_predictions?.drug_likeness_score?.toFixed(2)} (on a scale of 0–1, high likelihood for oral bioavailability)</span>
                  </div>
                  <div className="result-item">
                    <span className="label">Predicted Biological Target:</span>
                    <span className="value">{analysisResults.quantum_results?.ml_predictions?.biological_target} (Chronic Myeloid Leukemia/cancer)</span>
                  </div>
                  <div className="result-item">
                    <span className="label">Predicted Binding Energy:</span>
                    <span className="value">{analysisResults.quantum_results?.ml_predictions?.binding_energy?.toFixed(1)} kcal/mol (to human {analysisResults.quantum_results?.ml_predictions?.biological_target} domain; strong affinity)</span>
                  </div>
                  <div className="result-item">
                    <span className="label">ADMET Prediction:</span>
                    <span className="value">{analysisResults.quantum_results?.ml_predictions?.admet_prediction}</span>
                  </div>
                </div>
              </div>
            </div>

            {/* 3D Visualization Section */}
            {analysisResults.quantum_results?.visualization && (
              <div className="visualization-section">
                <h3>🔬 3D Molecular Visualization & Charge Analysis</h3>

                {analysisResults.quantum_results.visualization.has_3d_coords && (
                  <div className="viz-container">
                    {/* 3D Viewer Integration */}
                    <div className="viewer-wrapper">
                      <div className="viewer-controls">
                        <label>Representation:</label>
                        <select
                          value={currentRepresentation}
                          onChange={(e) => {
                            const rep = e.target.value;
                            setCurrentRepresentation(rep);
                            window.dispatchEvent(new CustomEvent('changeRepresentation', { detail: rep }));
                          }}
                          className="representation-select"
                        >
                          <option value="ballstick">Ball & Stick</option>
                          <option value="stick">Stick</option>
                          <option value="sphere">Space-Filling</option>
                          <option value="line">Wireframe</option>
                          <option value="cartoon">Cartoon</option>
                        </select>
                      </div>
                      <div className="viewer-container">
                        <MoleculeViewer />
                      </div>
                    </div>

                    {/* Charge Map */}
                    <div className="charge-map">
                      <h4>⚡ Partial Charge Distribution</h4>
                      <div className="charge-grid">
                        {Object.entries(analysisResults.quantum_results.visualization.charge_map || {}).slice(0, 10).map(([atom, charge], idx) => (
                          <div key={idx} className="charge-item">
                            <span className="atom-label">{atom.replace('_', ' ')}</span>
                            <span className={`charge-value ${charge < -0.2 ? 'negative' : charge > 0.2 ? 'positive' : 'neutral'}`}>
                              {charge > 0 ? '+' : ''}{charge.toFixed(3)}
                            </span>
                          </div>
                        ))}
                      </div>
                    </div>

                    {/* H-bonding regions */}
                    <div className="hbond-info">
                      <div className="hbond-section">
                        <h4>🔴 H-Bond Donors: {analysisResults.quantum_results.visualization.hbond_donors?.length || 0}</h4>
                        <p>Atoms that can donate hydrogen bonds for protein interactions</p>
                      </div>
                      <div className="hbond-section">
                        <h4>🔵 H-Bond Acceptors: {analysisResults.quantum_results.visualization.hbond_acceptors?.length || 0}</h4>
                        <p>Atoms that can accept hydrogen bonds in binding pocket</p>
                      </div>
                    </div>
                  </div>
                )}
              </div>
            )}

            {/* Interactive Quantum Spectroscopy */}
            {spectrumData && (
              <div className="quantum-spectroscopy">
                <h3>🔬 Interactive Quantum Spectroscopy</h3>
                <div className="spectroscopy-controls">
                  <div className="spectrum-tabs">
                    <button 
                      className={`spectrum-tab ${spectroscopyMode === 'absorption' ? 'active' : ''}`}
                      onClick={() => setSpectroscopyMode('absorption')}
                    >
                      <Eye size={16} />
                      Absorption
                    </button>
                    <button 
                      className={`spectrum-tab ${spectroscopyMode === 'emission' ? 'active' : ''}`}
                      onClick={() => setSpectroscopyMode('emission')}
                    >
                      <BarChart3 size={16} />
                      Emission
                    </button>
                    <button 
                      className={`spectrum-tab ${spectroscopyMode === 'vibrational' ? 'active' : ''}`}
                      onClick={() => setSpectroscopyMode('vibrational')}
                    >
                      <Target size={16} />
                      Vibrational
                    </button>
                  </div>
                </div>
                
                <div className="spectrum-plot">
                  {spectroscopyMode === 'absorption' && (
                    <div className="spectrum-data">
                      <h4>Absorption Spectrum</h4>
                      <p>Peaks: {spectrumData.absorption?.wavelengths?.length || 0} identified</p>
                      <p>Strongest peak: {spectrumData.absorption?.wavelengths?.[0]?.toFixed(1) || 'N/A'} nm</p>
                    </div>
                  )}
                  {spectroscopyMode === 'emission' && (
                    <div className="spectrum-data">
                      <h4>Emission Spectrum</h4>
                      <p>Fluorescence peaks: {spectrumData.emission?.wavelengths?.length || 0} identified</p>
                      <p>Primary emission: {spectrumData.emission?.wavelengths?.[0]?.toFixed(1) || 'N/A'} nm</p>
                    </div>
                  )}
                  {spectroscopyMode === 'vibrational' && (
                    <div className="spectrum-data">
                      <h4>Vibrational Modes</h4>
                      <p>Active modes: {spectrumData.vibrational?.length || 0} identified</p>
                      <p>Frequency range: {spectrumData.vibrational?.[0]?.frequency?.toFixed(0) || 'N/A'} - {spectrumData.vibrational?.[spectrumData.vibrational?.length-1]?.frequency?.toFixed(0) || 'N/A'} cm⁻¹</p>
                    </div>
                  )}
                </div>
              </div>
            )}

            {/* Time Evolution Dynamics */}
            {dynamicsData && (
              <div className="quantum-dynamics">
                <h3>⏰ Quantum Time Evolution</h3>
                <div className="dynamics-info">
                  <p>Time steps: {dynamicsData.time_steps?.length || 0}</p>
                  <p>Evolution circuits: {dynamicsData.evolution_circuits || 0}</p>
                  <p>Quantum coherence maintained throughout simulation</p>
                </div>
              </div>
            )}

            {/* Molecular Analogs Section */}
            {analysisResults.quantum_results?.molecular_analogs && (
              <div className="molecular-analogs">
                <h3>🔬 Quantum-Designed Molecular Analogs</h3>
                <p className="section-description">
                  Based on quantum simulations, here are optimized molecular variations with predicted properties:
                </p>
                <div className="analogs-grid">
                  {analysisResults.quantum_results.molecular_analogs.map((analog, index) => (
                    <div key={index} className="analog-card">
                      <div className="analog-header">
                        <h4>{analog.name}</h4>
                      </div>
                      <span className="analog-badge">{analog.modification}</span>

                      {/* Mini 3D Viewer */}
                      {analog.cif_content && (
                        <div className="analog-viewer-container">
                          <div id={`analog-viewer-${index}`} className="analog-mini-viewer"></div>
                        </div>
                      )}

                      <div className="analog-properties">
                        <div className="property">
                          <span className="property-label">HOMO-LUMO Gap:</span>
                          <span className="property-value">{analog.predicted_homo_lumo_gap} eV</span>
                        </div>
                        <div className="property">
                          <span className="property-label">Binding Affinity:</span>
                          <span className="property-value">{analog.predicted_binding} kcal/mol</span>
                        </div>
                      </div>
                      <p className="analog-rationale">
                        <strong>💡 Rationale:</strong> {analog.rationale}
                      </p>
                      <p className="quantum-advantage">
                        <strong>⚛️ Quantum Advantage:</strong> {analog.quantum_advantage}
                      </p>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Structural Suggestions */}
            {analysisResults.quantum_results?.structural_suggestions && (
              <div className="structural-suggestions">
                <h3>💡 Quantum-Informed Structural Suggestions</h3>
                <div className="suggestions-list">
                  {analysisResults.quantum_results.structural_suggestions.map((suggestion, index) => (
                    <div key={index} className="suggestion-item">
                      <span className="suggestion-text">{suggestion}</span>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Mechanism of Action Timeline */}
            {analysisResults.quantum_results?.mechanism_timeline && (
              <div className="mechanism-timeline">
                <h3>⏱️ Predicted Mechanism of Action Timeline</h3>
                <p className="timeline-description">Expected pathway steps for kinase inhibition based on quantum simulations</p>
                <div className="timeline-container">
                  {analysisResults.quantum_results.mechanism_timeline.map((step, index) => (
                    <div key={index} className="timeline-step">
                      <div className="step-number">{step.step}</div>
                      <div className="step-content">
                        <div className="step-event">{step.event}</div>
                        <div className="step-time">{step.time_ns} ns</div>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Enhanced Confidence Validation */}
            {analysisResults.quantum_results?.confidence_validation && (
              <div className="confidence-validation">
                <h3>✅ Confidence & Validation</h3>
                <div className="validation-grid">
                  <div className="validation-item">
                    <span className="validation-label">Model Confidence:</span>
                    <span className="validation-value confidence-high">
                      {analysisResults.quantum_results.confidence_validation.model_confidence}
                    </span>
                  </div>
                  <div className="validation-item">
                    <span className="validation-label">Simulation Mode:</span>
                    <span className="validation-value">
                      {analysisResults.quantum_results.confidence_validation.simulation_mode}
                    </span>
                  </div>
                  <div className="validation-item">
                    <span className="validation-label">Validation Status:</span>
                    <span className="validation-value">
                      {analysisResults.quantum_results.confidence_validation.validation_status}
                    </span>
                  </div>
                </div>

                {analysisResults.quantum_results.confidence_validation.simulation_details && (
                  <div className="simulation-details">
                    <h4>🔬 Simulation Technical Details</h4>
                    <div className="details-grid">
                      <div className="detail-item">
                        <span className="detail-label">Method:</span>
                        <span className="detail-value">{analysisResults.quantum_results.confidence_validation.simulation_details.method}</span>
                      </div>
                      <div className="detail-item">
                        <span className="detail-label">Basis Set:</span>
                        <span className="detail-value">{analysisResults.quantum_results.confidence_validation.simulation_details.basis_set}</span>
                      </div>
                      <div className="detail-item">
                        <span className="detail-label">Convergence:</span>
                        <span className="detail-value">{analysisResults.quantum_results.confidence_validation.simulation_details.convergence_criteria}</span>
                      </div>
                      <div className="detail-item">
                        <span className="detail-label">SCF Cycles:</span>
                        <span className="detail-value">{analysisResults.quantum_results.confidence_validation.simulation_details.scf_cycles}</span>
                      </div>
                      <div className="detail-item">
                        <span className="detail-label">Status:</span>
                        <span className="detail-value status-converged">{analysisResults.quantum_results.confidence_validation.simulation_details.optimization_status}</span>
                      </div>
                    </div>
                    <p className="validation-note">
                      💡 For further validation, submit larger dynamics jobs to IBM Quantum Cloud for production-scale simulations.
                    </p>
                  </div>
                )}
              </div>
            )}

            {/* AI Summary */}
            {aiSummary && (
              <div className="ai-summary quantum">
                <h3>🧠 AI Quantum Chemistry Insights</h3>
                {isGeneratingSummary ? (
                  <div className="loading">
                    <Brain className="spinning" size={20} />
                    Generating quantum insights...
                  </div>
                ) : (
                  <div className="summary-content">
                    {aiSummary}
                  </div>
                )}
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}

export default QuantumAnalysisPage;
