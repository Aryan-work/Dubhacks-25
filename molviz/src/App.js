import React, { useEffect } from 'react';
import { BrowserRouter, Routes, Route, Link } from 'react-router-dom';
import './App.css';
import Home from './pages/Home';
import ViewerPage from './pages/ViewerPage';
import QuantumAnalysisPage from './pages/QuantumAnalysisPage';
import statsigService from './services/statsigService';

function App() {
  useEffect(() => {
    // Initialize Statsig with session replay and auto-capture
    const initializeStatsig = async () => {
      try {
        await statsigService.initialize();
        statsigService.startSession();
        
        // Log app initialization
        statsigService.logEvent('app_initialized', {
          timestamp: new Date().toISOString(),
          version: '1.0.0',
          features: ['quantum_analysis', 'molecular_visualization', 'ai_insights']
        });
      } catch (error) {
        console.error('Failed to initialize Statsig:', error);
      }
    };

    initializeStatsig();

    // Cleanup on unmount
    return () => {
      statsigService.endSession();
    };
  }, []);

  const handleNavigation = (page) => {
    statsigService.logUserInteraction('navigation', 'navbar', { page });
  };

  return (
    <BrowserRouter>
      <div className="App">
        <nav className="nav">
          <Link to="/" onClick={() => handleNavigation('home')}>Home</Link>
          <Link to="/viewer" onClick={() => handleNavigation('viewer')}>Molecule Viewer</Link>
          <Link to="/quantum" onClick={() => handleNavigation('quantum')}>Quantum Chemistry</Link>
        </nav>
        <main className="main">
          <Routes>
            <Route path="/" element={<Home />} />
            <Route path="/viewer" element={<ViewerPage />} />
            <Route path="/quantum" element={<QuantumAnalysisPage />} />
          </Routes>
        </main>
      </div>
    </BrowserRouter>
  );
}

export default App;
