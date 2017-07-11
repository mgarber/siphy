package broad.pda.chromosome;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;


import broad.core.alignment.RepeatMaskerReader.RepeatStatistic;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GFF;
import broad.core.annotation.GenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.sequence.SequenceRegion;

public abstract class Mammal {
	protected File sequenceDirectory;
	private Chromosome X;
	private Chromosome Y;
	private Chromosome M;
	private List<Chromosome> Un;
	private ArrayList<Chromosome> autosomes;
	private ArrayList<Chromosome> random;
	private ArrayList<Chromosome> altHaplotypes;
	private File repeatDirectory;
	private File annotationDirectory;
	private Collection<RepeatStatistic> repeatStatistics;
	private Collection<RepeatStatistic> repeatFamilyStatistics;
	private Collection<RepeatStatistic> unpRepeatStatistics;
	private Collection<RepeatStatistic> unpRepeatFamilyStatistics;
	
	public abstract String getLetterCode();
	
	public Mammal(File sequenceDirectory) throws Exception {
		if(!sequenceDirectory.isDirectory()) {
			throw new Exception("sequence directory " + sequenceDirectory +" is not a directory");
		}
		this.sequenceDirectory = sequenceDirectory;
		String seqDirPath = sequenceDirectory.getAbsolutePath();
		seqDirPath = seqDirPath.lastIndexOf("/") == seqDirPath.length() -1 
			? seqDirPath.substring(0, seqDirPath.length() -1)
			: seqDirPath;
			
		this.repeatDirectory   = new File(seqDirPath + "_repeatinfo/");
		this.annotationDirectory = new File(seqDirPath + "_annotations/");
		this.autosomes = new ArrayList<Chromosome>();
		this.random = new ArrayList<Chromosome>();
		this.Un = new ArrayList<Chromosome>();
		this.altHaplotypes = new ArrayList<Chromosome>();
		
		String [] seqSubDirs = sequenceDirectory.list();		
		for(int i = 0; i < seqSubDirs.length; i++) {
			File dir = new File(sequenceDirectory.getAbsoluteFile() + "/" + seqSubDirs[i]);
			if(! dir.isDirectory() || dir.getName().contains("tmp")) {
				continue;
			}
			
			if(dir.getName().matches("Un[0-9]+")) {
				continue;
			}
			if(useMergedUn() && dir.getName().startsWith("Un") ) {
				continue;
			}
			//System.out.println("Processing " + dir + " looking for agp " + seqDirPath + "/" + seqSubDirs[i] + "/chr" + seqSubDirs[i] + ".agp");
			Chromosome chr = new Chromosome(seqDirPath + "/" + seqSubDirs[i] + "/chr" + seqSubDirs[i] + ".agp");
			if(chr.getSymbol().endsWith("_random") ) {
				random.add((chr));
			} else if("X".equals(seqSubDirs[i])) {
				//System.out.println("Setting X" );
				X = chr;
			} else if ("Y".equals(seqSubDirs[i])){
				//System.out.println("Setting Y");
				Y = chr;
			} else if ("M".equals(seqSubDirs[i])){
				//System.out.println("Setting M");
				M = chr;
			} else if (seqSubDirs[i].contains("Un")) {
				Un.add(chr);
				random.add(chr);
			} else if(seqSubDirs[i].contains("hap")){
				altHaplotypes.add(chr);
			}else {
				//System.out.println("Adding chromosome " + chr.getSymbol() + " to the autosomes");
				autosomes.add(chr);
			}
		}

		Collections.sort(autosomes, new Comparator<Chromosome>() {

			public int compare(Chromosome arg0, Chromosome arg1) {
				return (int) (arg1.getSize() - arg0.getSize());
			}
			
		});
	}
	
	public List<Chromosome> getAllNonRandomChromosomes() {
		ArrayList<Chromosome> all = new ArrayList<Chromosome>(autosomes.size() + 3);		
		all.addAll(autosomes);
		if(X != null) {
			all.add(X);
		}
		
		if(Y != null) {
			all.add(Y);
		}
		
		if(M != null) {
			all.add(M);
		}
		
		return all;
	}
	
	public List<Chromosome> getAllChromosomes() {
		List<Chromosome> all = getAllNonRandomChromosomes();
		all.addAll(random);
		
		return all;
	}
	
	public Chromosome getX() { return X;}
	
	public void loadChromosomeBands() throws IOException, ParseException {
		File bandAnnotation = new File(annotationDirectory + "/chromosome_bands.txt");
		if (bandAnnotation.exists()) {
			HashMap<String, List<GFF>> chrBandMap = new HashMap<String, List<GFF>>();
			String line;
			BufferedReader br = new BufferedReader(new FileReader(bandAnnotation));
			try {
				while((line = br.readLine()) != null) {
					if(line.startsWith("#")){
						continue;
					}
					//System.out.println(line.replace("\t", "-TAB-"));
					String [] lineSplit = line.split("\t");
					GFF gff = new GFF(lineSplit[3]);
					gff.setStart(Integer.parseInt(lineSplit[1]));
					gff.setEnd(Integer.parseInt(lineSplit[2]));
					gff.setChromosome(lineSplit[0].substring(3));
					gff.addAttribute("type", lineSplit[4]);

					List<GFF> sequenceGFFList = chrBandMap.get(gff.getChromosome());
					if(sequenceGFFList == null) {
						sequenceGFFList = new ArrayList<GFF>();
						chrBandMap.put(gff.getChromosome(), sequenceGFFList);
					}
					sequenceGFFList.add(gff);
				}
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				try {
					System.out.print("Closing "+bandAnnotation );
					br.close();
					System.out.print(" ..... Closed\n");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			Iterator<Chromosome> it = getAllNonRandomChromosomes().iterator();
			while(it.hasNext()) {
				Chromosome c = it.next();
				List<GFF> chrBands = chrBandMap.get(c.getSymbol());
				//System.out.println("Setting banding chr" + c.getSymbol()  + " had " + (chrBands != null ? chrBands.size() : "-") + " bands");
				c.setBanding(chrBands);
			}
		}
	}
	
	public void writeSizes(String outFile) throws IOException {
		BufferedWriter bw = new  BufferedWriter(new FileWriter(outFile));
		Iterator<Chromosome> chrIt = getAllChromosomes().iterator();
		while(chrIt.hasNext()) {
			Chromosome c = chrIt.next();
			bw.write("chr" + c.getSymbol() + "\t" + c.getSize());
			bw.newLine();
		}
		bw.close();
	}
	
	public long getNonRandomTotalSize() {
		long size = 0;
		Iterator<Chromosome> chrIt = getAllNonRandomChromosomes().iterator();
		while(chrIt.hasNext()) {
			Chromosome c = chrIt.next();
			System.out.println("Adding size of chr" + c.getSymbol() + ": " + c.getSize());
			size += c.getSize();
		}		
		return size;
	}
	
	public List<Chromosome> getAutosomes() {
		ArrayList<Chromosome> autosomesCopy = new ArrayList<Chromosome>(autosomes.size());
		autosomesCopy.addAll(autosomes);
		return autosomesCopy;
	}


	/**
	 * @param repeatStatMap
	 * @param tmpIt
	 */
	private void updateRepeatStatMap(HashMap<String, RepeatStatistic> repeatStatMap, Iterator<RepeatStatistic> statIt) {
		while(statIt.hasNext()) {
			RepeatStatistic chrStat = statIt.next();
			RepeatStatistic stat    = repeatStatMap.get(chrStat.getName());
			if(stat == null) {
				repeatStatMap.put(chrStat.getName(), chrStat);
			} else {
				stat.add(chrStat);
			}
		}
	}

	public int getAutosomeNumber() {
		return autosomes.size();
	}
	
	public Chromosome getChromosome(String chr) {
		String chrSymbol = chr;
		if(chrSymbol != null) {
			chrSymbol = chrSymbol.replace("chr", "");
		}
		Chromosome c = null;
		if(X != null && X.getSymbol().equals(chrSymbol)) {
			c = X;
		}else if((Y != null) && Y.getSymbol().equals(chrSymbol)) {
			c = Y;
		}else if((M != null) && M.getSymbol().equals(chrSymbol)) {
			c = M;
		} 
		
		if(c == null) {
			c = findChromosomeInList(autosomes, chrSymbol);
		} 
		
		if (c == null) {
			c = findChromosomeInList(random, chrSymbol);
		}
		return c;
	}
	
	public List<Chromosome> getUn() {
		return Un;
	}
	
	public long getGenomeSize() {
		long size = 0;
		if(X != null) {
			size += X.getSize();
		}
		
		if(Y != null) {
			size += Y.getSize();
		}
		
		Iterator<Chromosome> it = autosomes.iterator();
		while(it.hasNext()) {
			size += it.next().getSize();
		}
		return size;
	}
	
	public long getAutosomeGenomeSize() {
		long size = 0;

		Iterator<Chromosome> it = autosomes.iterator();
		while(it.hasNext()) {
			size += it.next().getSize();
		}
		return size;
	}
	
	/**
	 * Uses pseudo random number generation to draw a (non random) chromosome with
	 * chance proportional to the chromosome size
	 * @return Chromosome draw by random.
	 */
	public Chromosome  drawChromosome() {
		Chromosome chosen = null;
		double draw = Math.random();
		long genomeSize = getGenomeSize();
		long currentSize = 0;
		Iterator<Chromosome> it = getAllNonRandomChromosomes().iterator();
		while(it.hasNext()) {
			chosen = it.next();
			currentSize += chosen.getSize();
			if(currentSize/(double)genomeSize >= draw) {
				break;
			}
		}

		return chosen;
	}
	
	public SequenceRegion drawRandomRegion(int size) {
		Chromosome c = drawChromosome();
		return c.drawRandomRegion(size);
	}
	
	/**
	 * Uses pseudo random number generation to draw a (non random) autosome with
	 * chance proportional to the autosome size
	 * @return Chromosome draw by random.
	 */
	public Chromosome  drawAutosome() {
		Chromosome chosen = null;
		double draw = Math.random();
		long genomeSize = getAutosomeGenomeSize();
		long currentSize = 0;
		Iterator<Chromosome> it = getAutosomes().iterator();
		while(it.hasNext()) {
			chosen = it.next();
			currentSize += chosen.getSize();
			if(currentSize/(double)genomeSize >= draw) {
				break;
			}
		}

		return chosen;
	}
	
	public SequenceRegion drawRandomAutosomeRegion(int size) {
		Chromosome c = drawAutosome();
		return c.drawRandomRegion(size);
	}
	
	public void insertRegionAtRandom(SequenceRegion region) {
		Chromosome c = drawChromosome();
		c.insertAtRandomPoint(region);
	}
	
	public abstract File getPreferredSequenceDir();
	
	public File getSequenceDir() {
		return sequenceDirectory == null ? getPreferredSequenceDir() : sequenceDirectory;
	}
	
	public abstract String getName();
	/**
	 * to be overrided if repeats are in directories
	 * @return
	 */
	protected boolean areRepeatsInDeepDirHierarchy() {
		return !(new File(repeatDirectory.getAbsolutePath() + "/chr1.fa.out").exists());
	}
	
	protected boolean useMergedUn() {
		return false;
	}
	
	/**
	 * subclasses should override if the repeats contain overlapping segments (i.e. are merged with trf)
	 * @return
	 */
	protected boolean doRepeatListNeedsCleanup() {return true;}
	
	public File getRepeatFile(String chrSymbol) {
		return new File(repeatDirectory.getAbsolutePath() + "/"  + 
				(areRepeatsInDeepDirHierarchy() ? chrSymbol + "/" : "") + 
				"chr"+chrSymbol + ".fa.out");
	}
	
	public void loadRepeats(boolean loadStatisticsOnly) throws FileNotFoundException {
		Chromosome c = null;
		HashMap<String, RepeatStatistic> repeatStatMap = new HashMap<String, RepeatStatistic>();
		HashMap<String, RepeatStatistic> repeatFamStatMap = new HashMap<String, RepeatStatistic>();
		HashMap<String, RepeatStatistic> unpRepeatStatMap = new HashMap<String, RepeatStatistic>();
		HashMap<String, RepeatStatistic> unpRepeatFamStatMap = new HashMap<String, RepeatStatistic>();

		ArrayList<Chromosome> placedList = new ArrayList<Chromosome>(autosomes);
		if(X != null) {
			placedList.add(X);
		}
		Iterator<Chromosome> it = placedList.iterator();
		while(it.hasNext()) {
			c = it.next();
			File repeatFile = getRepeatFile(c.getSymbol());
			c.computeRepeatStatistics(repeatFile, loadStatisticsOnly);				
			updateRepeatStatMap(repeatStatMap, c.getRepeatStats().iterator());
			updateRepeatStatMap(repeatFamStatMap, c.getRepeatFamilyStats().iterator());
			if(loadStatisticsOnly) {
				c.unloadRepeats();
			}
		}
		repeatStatistics = repeatStatMap.values();
		repeatFamilyStatistics = repeatFamStatMap.values();
		
		it = random.iterator();
		while(it.hasNext()) {
			c = it.next();
			File repeatFile = new File(repeatDirectory.getAbsolutePath() + "/"  + 
					(areRepeatsInDeepDirHierarchy() ? c.getSymbol() + "/" : "") + 
					"chr"+c.getSymbol() + ".fa.out");
			if(repeatFile.exists()) {
				c.computeRepeatStatistics(repeatFile, doRepeatListNeedsCleanup());				
				updateRepeatStatMap(unpRepeatStatMap, c.getRepeatStats().iterator());
				updateRepeatStatMap(unpRepeatFamStatMap, c.getRepeatFamilyStats().iterator());
			}
			if(loadStatisticsOnly) {
				c.unloadRepeats();
			}
		}		
		unpRepeatStatistics = unpRepeatStatMap.values();
		unpRepeatFamilyStatistics = unpRepeatFamStatMap.values();
	} 
	
	private Chromosome findChromosomeInList(List<Chromosome> chrs, String symbol) {
		Iterator<Chromosome> it = chrs.iterator();
		Chromosome c = null;
		while(c == null && it.hasNext()) {
			Chromosome chr = it.next();
			if(symbol.equals(chr.getSymbol())) {
				c = chr;
			}
		}
		return c;
	}
	
	public Collection<RepeatStatistic> getRepeatStatistics() { return repeatStatistics;}
	public Collection<RepeatStatistic> getRepeatFamilyStatistics() { return repeatFamilyStatistics;}
	public Collection<RepeatStatistic> getUnplacedRepeatStatistics() { return unpRepeatStatistics;}
	public Collection<RepeatStatistic> getUnplacedRepeatFamilyStatistics() { return unpRepeatFamilyStatistics;}

	public File getAnnotationDirectory() {
		return annotationDirectory;
	}
	
	public List<GenomicAnnotation> getReducedRepresentation(int totalSize, List<GenomicAnnotation> toAvoid) {
		long totalGenomicSize = getGenomeSize();
		List<Chromosome> nonRandomChrs = getAutosomes();
		Iterator<Chromosome>  nonRandomChrIt = nonRandomChrs.iterator();
		ArrayList<GenomicAnnotation> chunks = new ArrayList<GenomicAnnotation>(nonRandomChrs.size());
		Chromosome c = null;
		while(nonRandomChrIt.hasNext()) {
			c = nonRandomChrIt.next();
			double chrSizePrct = c.getSize()/(double)totalGenomicSize;
			int chunkSize = (int) Math.round(totalSize * chrSizePrct);
			BasicGenomicAnnotation chunk = new BasicGenomicAnnotation("chr" + c.getSymbol() + "_chunk");
			chunk.setChromosome(c.getSymbol());
			chunk.setStart(1);
			chunk.setEnd(1 + chunkSize);
			chunks.addAll(c.shuffle(chunk, toAvoid));
			System.out.println("Added chunk of size " + chunk.getLength() + " for chr" + c.getSymbol() + " out of intended total " + totalSize);
		}
		
		return chunks;
	}

	public void delete(SequenceRegion region) {
		Chromosome c = getChromosome(region.getChromosome());
		c.delete(region);
	}

	public void invert(SequenceRegion region) {
		Chromosome c = getChromosome(region.getChromosome());
		c.invert(region);
	}
	
	public void insert(SequenceRegion region, GenomicAnnotation insertionPoint) {
		Chromosome c = getChromosome(insertionPoint.getChromosome());
		c.insert(region, insertionPoint);
	}

	public void translocate(SequenceRegion region1, SequenceRegion region2) {
		Chromosome c1 = getChromosome(region1.getChromosome());
		Chromosome c2 = getChromosome(region2.getChromosome());
		
		c1.delete(region1);
		GenomicAnnotation insertPoint1 = new BasicGenomicAnnotation(region1.getName(), region1.getChromosome(), region1.getStart(), region1.getStart() + 1);
		c1.insert(region2, insertPoint1);
		
		c2.delete(region2);
		GenomicAnnotation insertPoint2 = new BasicGenomicAnnotation(region2.getName(), region2.getChromosome(), region2.getStart(), region2.getStart() + 1);
		c2.insert(region1, insertPoint2);
	}

}
